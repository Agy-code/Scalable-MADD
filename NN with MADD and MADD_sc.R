# -------------------------------------------------------------------------
# Nearest Neighbor Classifier with MADD
# -------------------------------------------------------------------------
#
# This function implements NN_MADD and NN_MADD_sc classification methods.
#
# Steps:
#   1. For a query point z, compute rho_sc(z, x_i) for each training
#      observation x_i using a representative set X* (ref.points).
#      If X* is the full training set, this corresponds to NN_MADD.
#
#   2. For each class j, determine the minimum rho_sc(z, x_i)
#      among all training points belonging to class j.
#
#   3. Assign z to the class with the smallest class-wise minimum distance.
#
# Notes:
#   1. The function can be run in either memory-efficient or
#      computation-efficient mode:
#        - memory.eff = "YES":
#            Does not store the full n x n Euclidean distance matrix
#            (more memory-efficient, slightly slower).
#        - memory.eff = "NO":
#            Uses a precomputed or newly computed full distance matrix
#            (faster, but requires more memory).
#
#   2. print.MADD = TRUE returns the MADD matrix used in classification.
#      Default is FALSE.
#
# -------------------------------------------------------------------------
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


  class.labels <- unique(trainY)
  have_test_labels <- !is.null(testY)

  # -----------------------------------------------------------------------
  # Euclidean distances from every row of A to a single vector x
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

      # Use the given full distance matrix if available; otherwise compute it.
      if (is.null(Dtrain)) {
        Dtrain <- as.matrix(Rfast::Dist(trainX))
      } else {
        Dtrain <- as.matrix(Dtrain)
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

        # NN-classification:
        classwise_min <- sapply(class.labels, function(cl) {
          min(madd_row[trainY == cl])
        })

        preds[t] <- class.labels[which.min(classwise_min)]

        if (print.MADD) {
          MADD_mat[t, ] <- madd_row
        }
      }

    } else if (memory.eff == "YES") {

      # Memory-efficient branch:avoid storing the
      # full n x n training distance matrix.

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
  # This is the scalable version of MADD.
  # -----------------------------------------------------------------------
  } else {

    ref.points <- as.integer(ref.points)

    ref.points <- unique(ref.points)
    P <- length(ref.points)

    if (P < 1L) {
      stop("ref.points must contain at least one index.")
    }

    refX <- trainX[ref.points, , drop = FALSE]

    # Distances from every training point to every reference point.
    # If Dtrain is available, we reuse it; otherwise we compute only the
    # n x P block that we need here.
 
    if (!is.null(Dtrain)) {
      Dtrain <- as.matrix(Dtrain)
      D_tr_ref <- Dtrain[, ref.points, drop = FALSE]
    } else {
      D_tr_ref <- matrix(0, nrow = n, ncol = P)
      for (p in seq_len(P)) {
        D_tr_ref[, p] <- distance_to_point(trainX, refX[p, ])
      }
    }

    # Identify the training observations contained in the reference set.
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

  output <- list(
    MADD = if (print.MADD) MADD_mat else NULL,
    prediction = preds,
    conf.matrix = conf,
    accuracy = acc
  )

  return(output)
}
