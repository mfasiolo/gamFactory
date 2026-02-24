#' 
#' Function used in build_family_nl to check the correctness of g1 and g2
#' 
#' @name check_g_derivatives_fd
#' @rdname check_g_derivatives_fd
#' @export check_g_derivatives_fd
#' @description
#' only suitable for "si_nexp" effect
#' 


# find "nested" eff in olp$eff
.get_nested_index <- function(olp) {
  which(sapply(olp$eff, function(e) "nested" %in% class(e)))[1]
}

# get upper triangle index pairs like (1,1),(1,2),...,(1,na),(2,2),...,(na,na)
.upper_tri_pairs <- function(na) {
  idx_j <- integer(0); idx_k <- integer(0)
  for (j in seq_len(na)) {
    for (k in j:na) {
      idx_j <- c(idx_j, j)
      idx_k <- c(idx_k, k)
    }
  }
  data.frame(j = idx_j, k = idx_k)
}

# generate readable names for upper-tri using parameter names
.upper_tri_names <- function(par_names){
  out <- character()
  m <- length(par_names)
  for(i in 1:m) for(j in i:m) out <- c(out, paste0(par_names[i], " : ", par_names[j]))
  out
}

# ------------------------- main checker -------------------------
check_g_derivatives_fd <- function(olp, h1 = 1e-6, h2 = 2e-4, verbose = TRUE, per_sample = FALSE) {
  i_nested <- .get_nested_index(olp)
  if (is.na(i_nested)) stop("didn't find nested eff in olp$eff")
  e <- olp$eff[[i_nested]]
  if (is.null(e$eval) || !is.function(e$eval)) stop("didn't find e$eval")
  if (is.null(e$param)) stop("didn't find e$param")
  if (is.null(e$na)) stop("didn't find e$na (number of alpha params)")

  theta0 <- e$param
  na     <- e$na
  if(!is.null(e$alpha_center)) alpha_center <- e$alpha_center else alpha_center <- NULL
  if (length(theta0) < na) stop("param length < na")
  
  n      <- length(e$f)
  p_beta <- length(e$param) - na
  
  alpha0 <- theta0[seq_len(na)]
  beta0  <- if (p_beta > 0) theta0[(na + 1):length(theta0)] else numeric(0)
  
  # set up parameter names
  par_names <- if (!is.null(colnames(e$store$g1)) && ncol(e$store$g1) == na) {
    colnames(e$store$g1)
  } else {
    paste0("alpha[", seq_len(na), "]")
  }
  
  X_si <- e$store$X_si
  X_nexp <- e$store$X_nexp
  
  eval_g_at <- function(alpha_vec,alpha_center, X_si, X_nexp) {
    q <- ncol(X_nexp)
    p <- ncol(X_si)
    a0    <- exp(alpha_vec[1])
    anexp <- alpha_vec[2:(q + 1)]
    asi   <- alpha_vec[(q + 2):(1+p+q)]
    s     <- deriv_si_nexp(X_si, X_nexp, c(anexp, asi), deriv = 0,alpha_center = alpha_center)$d0
    a0 * (s - mean(s))
  }
  
  # ----------------- Gradient（per sample），fd with h1 -----------------
  grad_fd <- matrix(NA_real_, nrow = n, ncol = na)
  for (j in seq_len(na)) {
    ej <- rep(0, na); ej[j] <- 1
    g_plus  <- eval_g_at(alpha0 + h1 * ej, alpha_center, X_si, X_nexp)
    g_minus <- eval_g_at(alpha0 - h1 * ej, alpha_center, X_si, X_nexp)
    grad_fd[, j] <- (g_plus - g_minus) / (2 * h1)
  }
  colnames(grad_fd) <- par_names
  
  if (is.null(e$store$g1)) stop("can't find store$g1")
  grad_store <- e$store$g1
  if (!is.matrix(grad_store) || ncol(grad_store) != na) {
    stop("dim not match: store$g1 and na")
  }
  colnames(grad_store) <- par_names
  
  abs_diff_grad <- abs(grad_fd - grad_store)
  sum_grad <- c(
    max_abs    = max(abs_diff_grad),
    mean_abs   = mean(abs_diff_grad),
    median_abs = median(abs_diff_grad)
  )
  
  # ----------------- Hessian（per sample，upper triangle），with h2 -----------------
  if (is.null(e$store$g2)) stop("can't find store$g2")
  pairs <- .upper_tri_pairs(na)
  ntri  <- nrow(pairs)
  hess_store <- e$store$g2
  if (!is.matrix(hess_store) || ncol(hess_store) != ntri) {
    stop("dim not match: hess_store (g2) and upper-tri(na)")
  }
  ut_names <- .upper_tri_names(par_names)
  colnames(hess_store) <- ut_names
  
  g0 <- eval_g_at(alpha0, alpha_center, X_si, X_nexp)  # center point
  hess_fd <- matrix(NA_real_, nrow = n, ncol = ntri)
  for (col in seq_len(ntri)) {
    j <- pairs$j[col]
    k <- pairs$k[col]
    ej <- rep(0, na); ej[j] <- 1
    ek <- rep(0, na); ek[k] <- 1
    if (j == k) {
      # diag
      g_p <- eval_g_at(alpha0 + h2 * ej, alpha_center, X_si, X_nexp)
      g_m <- eval_g_at(alpha0 - h2 * ej, alpha_center, X_si, X_nexp)
      hess_fd[, col] <- (g_p - 2 * g0 + g_m) / (h2^2)
    } else {
      # non-diag
      g_pp <- eval_g_at(alpha0 + h2 * ej + h2 * ek, alpha_center, X_si, X_nexp)
      g_pm <- eval_g_at(alpha0 + h2 * ej - h2 * ek, alpha_center, X_si, X_nexp)
      g_mp <- eval_g_at(alpha0 - h2 * ej + h2 * ek, alpha_center, X_si, X_nexp)
      g_mm <- eval_g_at(alpha0 - h2 * ej - h2 * ek, alpha_center, X_si, X_nexp)
      hess_fd[, col] <- (fpp <- g_pp) - (fpm <- g_pm) - (fmp <- g_mp) + (fmm <- g_mm)
      hess_fd[, col] <- hess_fd[, col] / (4 * h2^2)
    }
  }
  colnames(hess_fd) <- ut_names
  
  abs_diff_hess <- abs(hess_fd - hess_store)
  sum_hess <- c(
    max_abs    = max(abs_diff_hess),
    mean_abs   = mean(abs_diff_hess),
    median_abs = median(abs_diff_hess)
  )
  
  # ----------------- combine matrix（set all sample's derivative together） -----------------
  # # gradient: (n*na)×2
  # grad_fd_vec    <- as.vector(t(grad_fd))      
  # grad_store_vec <- as.vector(t(grad_store))
  # # arrange by sample
  # grad_fd_vec    <- as.vector( t( t(grad_fd) ) )    
  # grad_store_vec <- as.vector( t( t(grad_store) ) )
  # combine by loop in order to make it clear
  grad_fd_vec <- grad_store_vec <- NULL
  grad_rn <- character()
  for (t in 1:n) {
    grad_fd_vec    <- c(grad_fd_vec,    grad_fd[t, ])
    grad_store_vec <- c(grad_store_vec, grad_store[t, ])
    grad_rn        <- c(grad_rn, paste0("s", t, ":", par_names))
  }
  grad_mat <- cbind(fd = grad_fd_vec, store = grad_store_vec)
  rownames(grad_mat) <- grad_rn
  
  # hessian: (n*ntri)×2，row name: "s{t}: alpha[i] : alpha[j]"
  hess_fd_vec    <- NULL
  hess_store_vec <- NULL
  hess_rn        <- character()
  for (t in 1:n) {
    hess_fd_vec    <- c(hess_fd_vec,    hess_fd[t, ])
    hess_store_vec <- c(hess_store_vec, hess_store[t, ])
    hess_rn        <- c(hess_rn, paste0("s", t, ":", ut_names))
  }
  hess_mat <- cbind(fd = hess_fd_vec, store = hess_store_vec)
  rownames(hess_mat) <- hess_rn
  
  result <- list(
    index_eff     = i_nested,
    deriv_used    = e$deriv,
    grad_fd       = grad_fd,
    grad_store    = grad_store,
    diff_grad     = abs_diff_grad,
    summary_grad  = sum_grad,
    hess_fd       = hess_fd,
    hess_store    = hess_store,
    diff_hess     = abs_diff_hess,
    summary_hess  = sum_hess,
    pairs         = pairs,
    par_names     = par_names,
    grad_mat      = grad_mat,   # combined matrix（gradient）
    hess_mat      = hess_mat    # combined matrix（hessian triangle）
  )
  
  # mean comparison
  if (verbose) {
    cat("== g1（gradient, column means）==\n")
    print(cbind(mean_fd = colMeans(grad_fd), mean_store = colMeans(grad_store)))
    cat("\n== g2（Hessian upper-tri, column means）==\n")
    print(cbind(mean_fd = colMeans(hess_fd), mean_store = colMeans(hess_store)))
    cat("\n")
  }
  
  # print combined matrix
  if (per_sample) {
      cat("Gradient (per-sample): FD vs Analytic\n")
      print(grad_mat, digits = 8)
      cat("\nHessian upper-tri (per-sample): FD vs Analytic\n")
      print(hess_mat, digits = 8)
  }
  
  invisible(result)
}
