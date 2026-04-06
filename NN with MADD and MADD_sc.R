# -------------------------------------------------------------------------
# Nearest Neighbor Classifier with MADD
# -------------------------------------------------------------------------
#
# This function implements NN_MADD and NN_MADD_sc classification methods.
#
# Steps:
#   1. For a query point z, compute rho_sc(z,.) for each training
#      observations using a representative set X* (ref.points).
#      If X* = NULL, then it uses the entire training data for MADD 
#      computations. Hence, it gives rho(z,.).
#
#   2. Then it applies the nearest neighbor classification technique 
#      based on the considered dissimilarity (rho or rho_sc). 
#
# Note: The function can be run in either memory-efficient or
#       computation-efficient mode:
#        - memory.eff = "YES":
#            Does not store the full n x n Euclidean distance matrix
#            (more memory-efficient, slightly slower).
#        - memory.eff = "NO":
#            Uses a full distance matrix
#            (faster, but requires more memory).
#
#  
# Input: 
#         trainX  = Features of training data, it should be a data.frame or matrix.
#         trainY  = Class labels of training data, it should be a vector.
#         testX   = Features of test data, it should be a data.frame or matrix.
#         testY   = Vector of class labels of test data, if available.
#         Dtrain  = Distance matrix between training data points. 
#                  This should be a distance object or a matrix.
#                   If not provided, uses the Rfast package.
#         ref.points = Reference points used in NN_MADD_sc classifier. 
#                      (There should be at least two reference points.) 
#         memory.eff = Memory efficiency needed or not. Default is NO.
#         print.MADD = Prints the MADD matrix in the output; it takes two 
#                      values, TRUE/FALSE. 
#                     - Default is FALSE.
# 
# Output:
#     -  MADD matrix, if print.MADD = TRUE.
#     - Predicted class-labels for the test data points.
#     - If true class labels of test data are available:
#         -- Confusion matrix.
#         -- Accuracy of the classifier.
# -------------------------------------------------------------------------

nn_madd <- function(trainX, 
                    trainY,
                    testX,
                    testY = NULL,
                    Dtrain = NULL,
                    ref.points = NULL,
                    memory.eff = "NO",
                    print.MADD = FALSE) {

  trainX <- as.matrix(trainX)
  testX  <- as.matrix(testX)
  trainY <- as.vector(trainY)

  if (!is.null(testY)) {
    testY <- as.vector(testY)
  }

  if (!memory.eff %in% c("NO", "YES")) {
    stop("memory.eff must be either 'NO' or 'YES'")
  }

  n <- nrow(trainX)    # Training data size
  m <- nrow(testX)     # Test data size
  d <- ncol(trainX)    # Dimension of the data.

  class.labels <- unique(trainY)
  have_test_labels <- !is.null(testY)

  # -----------------------------------------------------------------------
  # Euclidean distances from every row of A to a single vector x
  # -----------------------------------------------------------------------
  distance_to_point <- function(A, x) {
    sqrt(rowSums((A - matrix(x, nrow = nrow(A), ncol = ncol(A), byrow = TRUE))^2))
  }

  # -----------------------------------------------------------------------
  # Case 1: No reference point is supplied (calculates rho).
  # Then X* = trainX, i.e., the full training sample is used.
  # -----------------------------------------------------------------------
  if (is.null(ref.points)) {

    if (memory.eff == "NO") {

  # Uses the given full distance matrix if available; otherwise, computes it.
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
        d_z <- distance_to_point(trainX, z)

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

 # Memory-efficient branch: avoids storing the n x n distance matrix.

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
  # Case 2: Reference points are supplied (calculates rho_sc).
  # Then X* is the representative set indexed by ref.points.
  # -----------------------------------------------------------------------
  } else {
    P <- length(ref.points)
    
    refX <- trainX[ref.points, , drop = FALSE]
 
    if (!is.null(Dtrain)) {
      Dtrain <- as.matrix(Dtrain)
      D_tr_ref <- Dtrain[, ref.points, drop = FALSE]
    } else {
      D_tr_ref <- matrix(0, nrow = n, ncol = P)
      for (p in seq_len(P)) {
        D_tr_ref[, p] <- distance_to_point(trainX, refX[p, ])
      }
    }

    # To identify the training observations contained in the reference set.
    ref.position <- match(seq_len(n), ref.points, nomatch = 0L)
    is.reference <- ref.position > 0L

    if (print.MADD) {
      MADD_mat <- matrix(0, nrow = m, ncol = n)
    }

    preds <- character(m)

    for (t in seq_len(m)) {
      z <- testX[t, ]
      d_z_ref <- distance_to_point(refX, z)
      diff_mat <- abs(sweep(D_tr_ref, 2, d_z_ref, "-"))
      madd_row <- rowMeans(diff_mat)

      # If x_i itself belongs to the reference set, then it should be excluded
    
      if (any(is.reference)) {
        idx <- which(is.reference)
        col_idx <- ref.position[idx]
        adjusted_sum <- rowSums(diff_mat[idx, , drop = FALSE]) -
                                diff_mat[cbind(idx, col_idx)]
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
