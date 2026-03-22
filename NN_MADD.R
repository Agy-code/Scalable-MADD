# -------------------------------------------------------------------------
# Scalable MADD nearest-neighbor classifier
#
# This function implements the classification rule described in the paper:
#   1. Compute rho_sc(z, x_i) for every training observation x_i
#      using a representative set X* (called ref.points here).
#   2. For each class j, compute the minimum rho_sc(z, x_i)
#      over all training points belonging to class j.
#   3. Assign z to the class with the smallest classwise minimum.
#
# Important:
# - If ref.points is NULL, the full training set is used as the reference set.
# - memory.eff = "YES" avoids storing the full n x n training distance matrix.
# - memory.eff = "NO" uses a precomputed or newly computed full distance matrix
#   when ref.points is NULL.
# -------------------------------------------------------------------------

nn_madd <- function(trainX, trainY,
                             testX, testY = NULL,
                             Dtrain = NULL,
                             ref.points = NULL,
                             memory.eff = c("NO", "YES"),
                             print.MADD = FALSE) {

  memory.eff <- match.arg(memory.eff)

  trainX <- as.matrix(trainX)
  testX  <- as.matrix(testX)
  trainY <- as.vector(trainY)

  if (!is.null(testY)) {
    testY <- as.vector(testY)
  }

  n <- nrow(trainX)
  m <- nrow(testX)
  d <- ncol(trainX)

  if (ncol(testX) != d) {
    stop("trainX and testX must have the same number of columns.")
  }

  if (length(trainY) != n) {
    stop("Length of trainY must match the number of rows in trainX.")
  }

  class.labels <- unique(trainY)
  have_test_labels <- !is.null(testY)

  # -----------------------------------------------------------------------
  # Helper: Euclidean distances from every row of A to a single vector x
  # -----------------------------------------------------------------------
  distance_to_point <- function(A, x) {
    sqrt(rowSums((A - matrix(x, nrow = nrow(A), ncol = ncol(A), byrow = TRUE))^2))
  }

  # -----------------------------------------------------------------------
  # Case 1: No reference points supplied
  # Then X* = X, i.e. the full training sample is used as the reference set.
  # This reduces to ordinary MADD with the full training data.
  # -----------------------------------------------------------------------
  if (is.null(ref.points)) {

    if (memory.eff == "NO") {

      # Use the supplied full distance matrix if available; otherwise compute it.
      if (is.null(Dtrain)) {
        Dtrain <- as.matrix(Rfast::Dist(trainX))
      } else {
        Dtrain <- as.matrix(Dtrain)
      }

      if (!all(dim(Dtrain) == c(n, n))) {
        stop("Dtrain must be an n x n distance matrix.")
      }

      if (print.MADD) {
        MADD_mat <- matrix(0, nrow = m, ncol = n)
      }

      preds <- character(m)

      for (t in seq_len(m)) {
        z <- testX[t, ]

        # Distances from z to all training points
        d_z <- distance_to_point(trainX, z)

        # For each training point x_i, compute
        # mean_r | ||z - x_r|| - ||x_i - x_r|| |
        # over r != i
        diff_mat <- abs(sweep(Dtrain, 2, d_z, "-"))
        madd_row <- (rowSums(diff_mat) - diag(diff_mat)) / (n - 1)

        # Paper's rule:
        # For each class, take the minimum MADD to that class.
        classwise_min <- sapply(class.labels, function(cl) {
          min(madd_row[trainY == cl])
        })

        preds[t] <- class.labels[which.min(classwise_min)]

        if (print.MADD) {
          MADD_mat[t, ] <- madd_row
        }
      }

    } else if (memory.eff == "YES") {

      # Memory-efficient branch:
      # avoid storing the full n x n training distance matrix
      if (print.MADD) {
        MADD_mat <- matrix(0, nrow = m, ncol = n)
      }

      preds <- character(m)

      for (t in seq_len(m)) {
        z <- testX[t, ]
        d_z <- distance_to_point(trainX, z)

        madd_row <- numeric(n)

        for (i in seq_len(n)) {
          d_i <- distance_to_point(trainX, trainX[i, ])
          madd_row[i] <- mean(abs(d_z[-i] - d_i[-i]))
        }

        classwise_min <- sapply(class.labels, function(cl) {
          min(madd_row[trainY == cl])
        })

        preds[t] <- class.labels[which.min(classwise_min)]

        if (print.MADD) {
          MADD_mat[t, ] <- madd_row
        }
      }
    }

  # -----------------------------------------------------------------------
  # Case 2: Reference points supplied
  # Then X* is the representative set indexed by ref.points.
  # This is the scalable MADD used in the paper.
  # -----------------------------------------------------------------------
  } else {

    ref.points <- as.integer(ref.points)

    if (any(is.na(ref.points)) || any(ref.points < 1L) || any(ref.points > n)) {
      stop("ref.points must contain valid row indices of trainX.")
    }

    ref.points <- unique(ref.points)
    P <- length(ref.points)

    if (P < 1L) {
      stop("ref.points must contain at least one index.")
    }

    refX <- trainX[ref.points, , drop = FALSE]

    # Distances from every training point to every reference point.
    # If Dtrain is available, we reuse it; otherwise we compute only the
    # n x P block that we actually need.
    if (!is.null(Dtrain)) {
      Dtrain <- as.matrix(Dtrain)

      if (!all(dim(Dtrain) == c(n, n))) {
        stop("Dtrain must be an n x n distance matrix.")
      }

      D_tr_ref <- Dtrain[, ref.points, drop = FALSE]
    } else {
      D_tr_ref <- matrix(0, nrow = n, ncol = P)
      for (p in seq_len(P)) {
        D_tr_ref[, p] <- distance_to_point(trainX, refX[p, ])
      }
    }

    # Identify which training observations are themselves in the reference set.
    ref.position <- match(seq_len(n), ref.points, nomatch = 0L)
    is.reference <- ref.position > 0L

    if (print.MADD) {
      MADD_mat <- matrix(0, nrow = m, ncol = n)
    }

    preds <- character(m)

    for (t in seq_len(m)) {
      z <- testX[t, ]

      # Distances from z to the reference points
      d_z_ref <- distance_to_point(refX, z)

      # For each training point x_i:
      # rho_sc(z, x_i) = average over ref points r in X* \ {x_i}
      # of | ||z - r|| - ||x_i - r|| |
      diff_mat <- abs(sweep(D_tr_ref, 2, d_z_ref, "-"))
      madd_row <- rowMeans(diff_mat)

      # If x_i itself belongs to the reference set, then it should be excluded
      # from the average, exactly as in rho_sc(u, v) with X* \ {u, v}.
      if (any(is.reference)) {
        idx <- which(is.reference)
        col_idx <- ref.position[idx]

        adjusted_sum <- rowSums(diff_mat[idx, , drop = FALSE]) -
          diff_mat[cbind(idx, col_idx)]

        if (P == 1L) {
          stop("When a training point is also a reference point, at least two reference points are needed.")
        }

        madd_row[idx] <- adjusted_sum / (P - 1)
      }

      classwise_min <- sapply(class.labels, function(cl) {
        min(madd_row[trainY == cl])
      })

      preds[t] <- class.labels[which.min(classwise_min)]

      if (print.MADD) {
        MADD_mat[t, ] <- madd_row
      }
    }
  }

  conf <- NULL
  acc  <- NULL

  if (have_test_labels) {
    conf <- table(prediction = preds, truth = testY)
    acc  <- sum(diag(conf)) / sum(conf)
  }

  out <- list(
    MADD = if (print.MADD) MADD_mat else NULL,
    prediction = preds,
    conf.matrix = conf,
    accuracy = acc
  )

  return(out)
}
