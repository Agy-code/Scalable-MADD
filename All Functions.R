#--------------------------------------------------------------------------------------------------------
#                    (1) MADD_knn classifier (supports prototype-based method)
#--------------------------------------------------------------------------------------------------------

predict_madd_knn <- function(trainX, trainY,
                             testX, testY   = NULL,
                             Dtrain         = NULL, 
                             k              = 1,         # k-nn number of neighbors
                             prototype.ind  = NULL,
                             n_tol          = 6000L,
                             print.MADD     = FALSE) {

  # Dtrain = Distance matrix of the training data points, 
  #          which can be provided; if NULL and n <= n_tol,
  #          it will be computed via Rfast::Dist(trainX).

  trainX <- as.matrix(trainX)
  trainY <- as.vector(trainY)
  testX  <- as.matrix(testX)

  n <- nrow(trainX) # Training sample size
  m <- nrow(testX)  # Test sample size
  d <- ncol(trainX) # Dimension

  have_test_labels <- !is.null(testY)
  if (have_test_labels) testY <- as.vector(testY)

  # --------------------------------------------------------
  # CASE 1: no prototypes, n <= n_tol  
  # --------------------------------------------------------

  if (is.null(prototype.ind) && n <= n_tol) {

    # If Dtrain not supplied

    if (is.null(Dtrain))  Dtrain <- Rfast::Dist(trainX)      
   
    Dtrain <- as.matrix(Dtrain)      # Full n x n matrix

    if (print.MADD ) {
      MADD_mat <- matrix(0, nrow = m, ncol = n)
    }

    preds <- character(m)

    for (t in seq_len(m)) {
      z <- testX[t, ]

      d_z <- sqrt(rowSums((trainX - matrix(z, nrow = n, ncol = d, byrow = TRUE))^2))
      diff_mat <- abs(sweep(Dtrain, 2, d_z, "-"))
      row_sums <- rowSums(diff_mat) - diag(diff_mat)
      madd_row <- row_sums / (n - 1)

      # k-NN for this test point
      nn <- order(madd_row)[seq_len(k)]
      preds[t] <- names(which.max(table(trainY[nn])))

      if (print.MADD ) {
        MADD_mat[t, ] <- madd_row
      }
    }

  # --------------------------------------------------------
  # CASE 2: no prototypes, n > n_tol  
  # --------------------------------------------------------
  } else if (is.null(prototype.ind) && n > n_tol) {

    if (print.MADD ) {
      MADD_mat <- matrix(0, nrow = m, ncol = n)
    }
    preds <- character(m)

    for (t in seq_len(m)) {
      z <- testX[t, ]
      z<-as.numeric(z)
      d_z <- sqrt(rowSums((trainX - matrix(z, nrow = n, ncol = d, byrow = TRUE))^2))

      madd_row <- numeric(n)
      for (i in seq_len(n)) {
        d_i <- sqrt(rowSums((trainX - matrix(trainX[i, ], nrow = n, ncol = d, byrow = TRUE))^2))
        madd_row[i] <- mean(abs(d_z[-i] - d_i[-i]))
      }

      nn <- order(madd_row)[seq_len(k)]
      preds[t] <- names(which.max(table(trainY[nn])))

      if (print.MADD ) {
        MADD_mat[t, ] <- madd_row
      }
    }

 # --------------------------------------------------------
  # CASE 3: prototype-based
  # --------------------------------------------------------
  } else {
    prototype.ind <- as.integer(prototype.ind)
    P <- length(prototype.ind)
    protoX <- trainX[prototype.ind, , drop = FALSE]

 if (!is.null(Dtrain)) {
      # Use provided full distance matrix
      D_full <- as.matrix(Dtrain)               # n x n
      D_tp   <- D_full[, prototype.ind, drop = FALSE]  # n x P
    } else {
      # Direct computation
      D_tp <- matrix(0, nrow = n, ncol = P)
      for (j in seq_len(P)) {
        pj <- protoX[j, ]
        D_tp[, j] <- sqrt(rowSums(
          (trainX - matrix(pj, nrow = n, ncol = d, byrow = TRUE))^2
        ))
      }
    }
    proto_pos <- match(seq_len(n), prototype.ind, nomatch = 0L)
    is_proto  <- proto_pos > 0L

    if (print.MADD ) {
      MADD_mat <- matrix(0, nrow = m, ncol = n)
    }
    preds <- character(m)

    for (t in seq_len(m)) {
      z <- testX[t, ]
      d_zp <- sqrt(rowSums((protoX - matrix(z, nrow = P, ncol = d, byrow = TRUE))^2))

      diff_mat <- abs(sweep(D_tp, 2, d_zp, "-"))  # n x P

      base_sum  <- rowSums(diff_mat)
      base_mean <- base_sum / P

      if (any(is_proto)) {
        idx_proto_rows <- which(is_proto)
        col_proto_rows <- proto_pos[is_proto]
        adj_sum  <- base_sum[is_proto] - diff_mat[cbind(idx_proto_rows, col_proto_rows)]
        adj_mean <- adj_sum / (P - 1)
        base_mean[is_proto] <- adj_mean
      }

      madd_row <- base_mean

      nn <- order(madd_row)[seq_len(k)]
      preds[t] <- names(which.max(table(trainY[nn])))

      if (print.MADD ) {
        MADD_mat[t, ] <- madd_row
      }
    }
  }


  # --------------------------------------------------------
  # confusion matrix + accuracy if labels given
  # --------------------------------------------------------
  conf <- NULL
  acc  <- NULL
  if (have_test_labels) {
    conf <- table(preds, testY)
    acc  <- sum(diag(conf)) / sum(conf)
  }

  if (!print.MADD ) {
    list(MADD = NULL, prediction = preds, 
         conf.matrix = conf, accuracy = acc)
  } else {
    list(MADD = MADD_mat, prediction = preds,
          conf.matrix = conf, accuracy = acc)
  }
}




#--------------------------------------------------------------------------------------------------------
#                                   (2) DPP Sampler
#--------------------------------------------------------------------------------------------------------


#---------------------(i) DPP sampling from given eigen decomposition------------------

sample_DPP <- function(eigendecomposition, 
                       method = c("standard", "kDPP"), 
                       k = NULL  #Required for kDPP.
                       ) 
{  
  method <- match.arg(method)
  if (is.null(k)&& method == "kDPP") stop("Please provide 'k' for method = 'kDPP'")
  vals <- eigendecomposition$values
  vecs <- eigendecomposition$vectors
  N <- length(vals)


  orthogonalize_to_e <- function(V, e, tol = 1e-10) {
  V <- as.matrix(V)
  e <- as.matrix(e)

  # --- Step 1. Find a column v0 with non-zero inner product with e ---
  v_0 <- NULL
  for (j in seq_len(ncol(V))) {
    inner <- drop(t(V[, j, drop = FALSE]) %*% e)
    if (abs(inner) > tol) {
      v_0 <- V[, j, drop = FALSE]
      break
    }
  }

  # --- If, all columns are (numerically) orthogonal to e ---
  # Then span(V) is the required subspace.

  if (is.null(v_0)) {
    return(qr.Q(qr(V)))
  }

  # --- Step 2. Build W: Vectors in V but orthogonal to e ---
  W <- apply(V, 2, function(v) {
    v <- as.matrix(v)
    cont <- drop(t(v) %*% e / (t(v_0) %*% e))
    v - cont * v_0
  })

  W <- as.matrix(W)

  # --- 4. Drop (numerically) zero columns ---
  nonzero_cols <- colSums(abs(W)) > tol
  W <- W[, nonzero_cols, drop = FALSE]

  # If everything vanished, the intersection is {0}, return 0-column matrix
  if (ncol(W) == 0L) {
    return(matrix(0, nrow = nrow(V), ncol = 0L))
  }

  # --- 5. Orthonormal basis for span(W)---
  d <- qr(W)
  Q <- qr.Q(d)
  orthogonal_basis <- Q[, seq_len(d$rank), drop = FALSE]
  return(orthogonal_basis)
}


compute_elementary_symmetric <- function(k, eigenvalues) {
    N <- length(eigenvalues)
    e <- matrix(0, nrow = N + 1, ncol = k + 1)
    e[, 1] <- 1
    for (l in 2:(k + 1)) {
      for (n in 2:(N + 1)) {
        e[n, l] <- e[n - 1, l] + eigenvalues[n - 1] * e[n - 1, l - 1]
      }
    }
    return(e)
  }


sample_k_eigenvectors <- function(k, eigenvalues) {
  N <- length(eigenvalues)

  if (k > N) stop("k must be <= length(eigenvalues)")
  if (k < 0) stop("k must be non-negative")

  e_k_n <- compute_elementary_symmetric(k, eigenvalues)

  J <- integer(0)
  l <- k + 1

  for (n in (N + 1):2) {
    if (l == 1L) break

    denom <- e_k_n[n, l]
    if (denom == 0) {
      prob <- 0
    } else {
      prob <- eigenvalues[n - 1L] * e_k_n[n - 1L, l - 1L] / denom
    }

    prob <- max(0, min(1, prob))

    if (runif(1) < prob) {
      J <- c(J, as.integer(n - 1L))
      l <- l - 1L
    }
  }
  return(J)
}


if (method == "standard") {
    p_vals <- vals / (vals + 1)
    J <- which(runif(length(vals)) < p_vals)
    V <- as.matrix(vecs[, J, drop = FALSE])
  } else{
    J <- sample_k_eigenvectors(k, vals)
    V <- as.matrix(vecs[, J, drop = FALSE])
  } 

  # Sampling set Y
  Y <- integer(0)
  while (ncol(V) > 0) {
    sample_prob <- rowSums(V^2) / ncol(V)
    i <- sample(1:nrow(V), 1, prob = sample_prob)
    Y <- c(Y, as.integer(i))
    e_i <- rep(0, nrow(V))
    e_i[i] <- 1
    if (ncol(V) == 1) {
      break
    } else {
      V <- orthogonalize_to_e(V, as.matrix(e_i))
    }
  }
 return(as.vector(Y))
}


#---------------------(ii) DPP sampler for classification problem------------------

# It uses RFF for approximating DPP when n is huge.


dpp_index <- function(trainY, Dtrain, 
                      sample_size   = function(n) ceiling(sqrt(n)),
                      k             = NULL,
                      sample_method = "kDPP",
                      trainX        = NULL, #trainX is needed for approximate DPP
                      D_RFF         = 500L,
                      n_tolerance   = 2000L,
                      output        = "vector") 
{
  classes <- unique(trainY)
  
  # Ensure every class has at least 2 samples
  if (any(table(trainY) <= 1L)) {
    stop("At least one class has <= 1 member!")
  }

  all_ids   <- integer(0)                        # for "vector" mode
  class_ids <- vector("list", length(classes))   # for "list" mode
  names(class_ids) <- classes

  # --- helper: RFF using cos+sin, 1/sqrt(D) scaling ---
  RFF <- function(X, D, sigma) {
    n <- nrow(X)
    d <- ncol(X)
    # w_l ~ N(0, I_d / sigma^2)
    W <- matrix(rnorm(D * d, sd = 1 / sigma), nrow = D)  # D x d
    proj <- X %*% t(W)        # n x D
    Z <- cbind(cos(proj), sin(proj)) / sqrt(D)  # n x (2D)
    Z
  }

  for (cl in classes) {
    idx_cl <- which(trainY == cl)
    n_cl   <- length(idx_cl)

    # how many to sample for this class?
    k_cl <- if (!is.null(k)) {
      min(k, n_cl)
    } else {
      min(sample_size(n_cl), n_cl)
    }

    # ------------------------
    # CASE A: exact DPP if n_cl <= n_tolerance
    # ------------------------
    if (n_cl <= n_tolerance) {
      if (is.null(Dtrain)) {
        stop("Dtrain must be provided for exact DPP mode (n_cl <= n_tolerance).")
      }
      Dcl   <- Dtrain[idx_cl, idx_cl, drop = FALSE]
      sigma <- median(Dcl[upper.tri(Dcl)])

      Lcl <- exp(-(Dcl^2) / (2 * sigma^2))

      eig <- eigen(Lcl, symmetric = TRUE)
      sampled_local <- sample_DPP(eig, method = sample_method, k = k_cl)

    } else {

      # ------------------------
      # CASE B: RFF-based dual DPP if n_cl > n_tolerance
      # ------------------------
      if (!is.null(Dtrain)) {
        # if Dtrain is available and not too huge, we can still use its subset
        Dcl <- Dtrain[idx_cl, idx_cl, drop = FALSE]
        # maybe subsample above some size, but simple version:
        sigma <- median(Dcl[upper.tri(Dcl)])
      } else {
        # fallback: use distances in feature space on a subset
        train_cl <- trainX[idx_cl, , drop = FALSE]
        samp <- sample(n_cl, min(2000L, n_cl))
        D_samp <- as.matrix(dist(train_cl[samp, , drop = FALSE]))
        sigma <- median(D_samp[upper.tri(D_samp)])
      }

      train_cl <- as.matrix(trainX[idx_cl, , drop = FALSE])

      Z        <- RFF(X = train_cl, D = D_RFF, sigma = sigma)  # n_cl x (2D_RFF)
      C <- crossprod(Z)   # (2D_RFF) x (2D_RFF)

      eigC  <- eigen(C, symmetric = TRUE)
      vals  <- eigC$values
      V_feat <- eigC$vectors      # (2D_RFF) x r

      U <- Z %*% V_feat           # n_cl x r
      # normalize columns of U:
      U <- sweep(U, 2, sqrt(colSums(U^2)), "/")

      eig <- list(values = vals, vectors = U)
      r_eff<-sum(eig$values>1e-10)
      if(k_cl>r_eff) print("Too many samples! Adjusting the sample size")
      sampled_local <- sample_DPP(eig, method = sample_method, k = min(k_cl,r_eff))
    }

    sampled_global <- idx_cl[sampled_local]

    all_ids <- c(all_ids, sampled_global)
    class_ids[[as.character(cl)]] <- sampled_global
  }
  
  if (output == "vector") {
    return(sort(unique(all_ids)))
  } else {
    return(class_ids)
  }
}

#----------------- To find effective number of samples that can be drawn------------------

Find_reff <- function(trainX, trainY,
                      Dtrain = NULL,
                      D_RFF  = 500L,
                      tol    = 1e-10) {
  
   classes <- sort(unique(trainY))
  
  
  # --- helper: RFF using cos+sin, 1/sqrt(D) scaling ---
  RFF <- function(X, D, sigma) {
    n <- nrow(X)
    d <- ncol(X)
    # w_l ~ N(0, I_d / sigma^2)
    W <- matrix(rnorm(D * d, sd = 1 / sigma), nrow = D)  # D x d
    proj <- X %*% t(W)                                  # n x D
    Z <- cbind(cos(proj), sin(proj)) / sqrt(D)          # n x (2D)
    Z
  }
  
  r_eff_vec <- numeric(length(classes))
  names(r_eff_vec) <- classes
  
  for (cl in classes) {
    idx_cl <- which(trainY == cl)
    n_cl   <- length(idx_cl)
    
    # Training data for this class
    train_cl <- as.matrix(trainX[idx_cl, , drop = FALSE])
    
    # Estimate sigma
    if (!is.null(Dtrain)) {
      Dcl   <- Dtrain[idx_cl, idx_cl, drop = FALSE]
      sigma <- median(Dcl[upper.tri(Dcl)])
    } else {
      samp  <- sample(n_cl, min(2000L, n_cl))
      D_samp <- as.matrix(dist(train_cl[samp, , drop = FALSE]))
      sigma  <- median(D_samp[upper.tri(D_samp)])
    }
    
      Z <- RFF(X = train_cl, D = D_RFF, sigma = sigma)  # n_cl x (2D_RFF)
      C <- crossprod(Z)   # (2D_RFF) x (2D_RFF)

      eigC  <- eigen(C, symmetric = TRUE)
      vals  <- eigC$values
      r_eff <- sum(vals>tol)
  
    r_eff_vec[as.character(cl)] <- r_eff
  }  
  return(r_eff_vec)
}



#------------------------------------------------------------------------------------------------------------
#                 (3) Cross validation to choose number of prototypes
#------------------------------------------------------------------------------------------------------------


cross_validation <- function(trainX, trainY, Dtrain, k_choices,
                             folds = 5, knn = 1,D_RFF = 500L,
                             method = "kDPP") {
  Dtrain    <- as.matrix(Dtrain)
  kmax_req  <- max(k_choices)

  # simple k-NN on a distance block (rows=test, cols=train)
  knn_predict <- function(dist_block, labels, k = 1L) {
    if (k == 1L) {
      labels[max.col(-dist_block, ties.method = "first")]
    } else {
      apply(dist_block, 1L, function(r) {
        idx  <- head(order(r), k)
        labs <- labels[idx]
        names(which.max(table(labs)))
      })
    }
  }

  cv_folds_id <- caret::createFolds(trainY, k = folds, returnTrain = FALSE)
  results <- list()

  for (f in seq_along(cv_folds_id)) {
    te <- cv_folds_id[[f]]
    tr <- setdiff(seq_len(nrow(trainX)), te)
    trainX_f<-  as.matrix(trainX)[tr, ]
    trainY_f <- trainY[tr]
    n_train  <- length(tr)
    n_test   <- length(te)

    # Distance blocks from full Dtrain
    Dtrain_f       <- Dtrain[tr, tr, drop = FALSE]  # (n_train x n_train)
    D_test_train_f <- Dtrain[te, tr, drop = FALSE]  # (n_test  x n_train)

    # Sample per-class prototypes up to kmax_req
    dpp_ids <- dpp_index(trainX= trainX_f, trainY = trainY_f, Dtrain = Dtrain_f,
                         k = kmax_req, sample_method = method, output = "list")

    kmax_fold <- min(sapply(dpp_ids, length))
    k_use     <- k_choices[k_choices <= kmax_fold]
    if (length(k_use) == 0L) next

    # ---- Incremental state ----
    S_sum       <- matrix(0.0, nrow = n_test, ncol = n_train)
    T_added     <- 0L
    is_included <- integer(n_train)

    add_prototypes <- function(new_ids) {
      if (length(new_ids) == 0L) return(invisible(NULL))
      m_new   <- length(new_ids)
      dzx_mat <- D_test_train_f[, new_ids, drop = FALSE]   # n_test × m_new
      dxx_mat <- Dtrain_f[, new_ids, drop = FALSE]         # n_train × m_new

      diff_sums <- vapply(seq_len(n_test), function(i) {
                           diffs <- abs(
                 matrix(rep(dzx_mat[i, ], n_train), byrow = TRUE, ncol = m_new) - dxx_mat
                                        )
                          rowSums(diffs)
                           }, FUN.VALUE = numeric(n_train))

      S_sum <<- S_sum + t(diff_sums)
      T_added <<- T_added + m_new

      for (jproto in new_ids) {
        S_sum[, jproto] <<- S_sum[, jproto] - abs(
          D_test_train_f[, jproto] - Dtrain_f[jproto, jproto]
        )
        is_included[jproto] <<- 1L
      }
    }

    evaluate_now <- function() {
      denom <- pmax(1L, T_added - is_included)
      MADD_mean <- S_sum / matrix(denom, nrow = n_test, ncol = n_train, byrow = TRUE)
      pred_mean <- knn_predict(MADD_mean, trainY_f, k = knn)
      mean(pred_mean == trainY[te])
    }

    # ---- First k ----
    k_prev    <- k_use[1]
    first_ids <- unlist(lapply(dpp_ids, function(ix) ix[seq_len(k_prev)]), use.names = FALSE)
    add_prototypes(first_ids)

    m <- evaluate_now()
    results[[length(results) + 1L]] <- data.frame(
      fold = f,
      k_per_class = k_prev,
      accuracy = m,
      stringsAsFactors = FALSE
    )

    # ---- Remaining k's ----
    if (length(k_use) > 1L) {
      for (k_cur in k_use[-1]) {
        new_ids <- unlist(lapply(dpp_ids, function(ix) ix[(k_prev + 1):k_cur]), use.names = FALSE)
        add_prototypes(new_ids)

        m <- evaluate_now()
        results[[length(results) + 1L]] <- data.frame(
          fold = f,
          k_per_class = k_cur,  
          accuracy = m,
          stringsAsFactors = FALSE
        )
        k_prev <- k_cur
      }
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#----------------------------------------------- x x x x x--------------------------------------------------------------
