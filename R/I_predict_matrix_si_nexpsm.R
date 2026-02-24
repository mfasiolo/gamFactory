#'
#' Predict using nested exponential smoothing effects
#' 
#'
#' @noRd

.predict.matrix.si_nexpsm <- function(object, data, get.xa = FALSE){
  
  # -----------------1, prepare information--------------------
  si <- object$xt$si
  alpha <- si$alpha
  positive_si <- isTRUE(si$positive_si)
  
  n_si <- si$n_si
  n_nexp <- si$n_nexp
  
  original_data <- data[[object$term]]
  times <- NULL
  nms  <- colnames(original_data)
  i_t  <- which(nms == "times")            # index of 'times' column, if any
  times <- NULL                            # extract time column if present
  if (length(i_t)) {
    times <- original_data[ , i_t]
    Xall  <- original_data[ , -i_t, drop = FALSE]
  }
  
  if(length(i_t)){
    all_columns <- ncol(Xall)
    all_columns_true <- n_si + n_nexp
    nrep <- ceiling(all_columns / all_columns_true)
    n_obs  <- nrow(Xall)
    stopifnot(ncol(Xall) == nrep * (n_si + n_nexp))
    
    X_block <- Xall[, 1:(nrep * n_si), drop = FALSE]
    W_block <- Xall[, (nrep * n_si + 1):(nrep * (n_si + n_nexp)), drop = FALSE]
    
    X_si_long   <- matrix(NA_real_, nrow = n_obs * nrep, ncol = n_si)
    X_nexp_long <- matrix(NA_real_, nrow = n_obs * nrep, ncol = n_nexp)
    
    # recover
    for (d in 1:n_obs) {
      for (s in 1:nrep) {
        row <- (d - 1) * nrep + s
        X_si_long[row, ]   <- X_block[d, ((s - 1) * n_si + 1):(s * n_si)]
        X_nexp_long[row, ] <- W_block[d, ((s - 1) * n_nexp + 1):(s * n_nexp)]
      }
    }
    
    # names（eg. c("East","South","West","North")、c("intercept","time_diff")）
    colnames(X_si_long)   <- paste0("x", 1:n_si)      #  c("East","South","West","North")
    colnames(X_nexp_long) <- paste0("w", 1:n_nexp)    #  c("intercept","time_diff")
    
    # --- original X_si, X_nexp---
    X_si   <- X_si_long
    X_nexp <- X_nexp_long
    
    # Drop INF
    ok_si   <- rowSums(is.finite(X_si))   == ncol(X_si)
    ok_nexp <- rowSums(is.finite(X_nexp)) == ncol(X_nexp)

    last_good <- max(which(ok_si & ok_nexp))

    X_si   <- X_si[  1:last_good, , drop = FALSE]
    X_nexp <- X_nexp[1:last_good, , drop = FALSE]
  }else{
    
    X_si <- original_data[, 1:n_si, drop = FALSE]
    X_nexp <- original_data[, (n_si+1):(n_si+n_nexp), drop = FALSE]
  }
  
  # alpha =[alpha_nexp, alpha_si]
  alpha_nexp <- alpha[1:n_nexp]  # without scale parameter
  alpha_si <- alpha[(n_nexp + 1):(n_nexp + n_si)]
  alpha_center <- si$alpha_center
  # sd_g <- si$sd_g
  
  
  # -----------------2, deal with xm and center X_si--------------------
  # original xm only contains mean of X_si, we have to store the mean of X_nexp
  # currently use 0
  if(length(si$xm) == n_si){
    si$xm <- c(si$xm, nexp = 0)
  }
  
  if (isTRUE(si$center)){
    X_si <- scale(X_si, center = si$xm[1:n_si], scale = FALSE)
  } 
  
  # -----------------3, rescaling by B matrix, get z--------------------
  X_si <- X_si %*% si$B_si
  X_nexp <- X_nexp %*% si$B_nexp
  
  res <- deriv_si_nexp(X_si = X_si,
                           X_nexp = X_nexp,
                           param = c(alpha_nexp, alpha_si),
                           times = times,
                           deriv = 1,
                           alpha_center = alpha_center,
                           Z0 = si$Z0,
                           positive_si = positive_si)
  xsm_raw <- res$d0
  xsm_unscaled <- xsm_raw - si$xm["nexp"]  # center
  xsm <- xsm_unscaled  # scaling
  # xsm <- xsm_unscaled/sd_g
  
  # ---------------4, first driv--------------------
  xa_da <- res$d1  #n × (n_si+n_nexp)
  # xa_da <- dres$d1/sd_g
  
  if(get.xa){
    # construct derivatives manually
    
    # ∂s/∂alpha_nexp
    xa_dalpha_nexp <- xa_da[, 1:n_nexp, drop = FALSE]  # n × n_nexp
    
    # ∂s/∂alpha_si
    xa_dalpha_si <- xa_da[, (n_nexp + 1):(n_si + n_nexp), drop = FALSE]  # n × n_si
    
    return(list(
      X = NULL,  # don't need it for now, because type = "link" will call basis again
      xa = xsm,
      xa_dalpha_si = xa_dalpha_si,
      xa_dalpha_nexp = xa_dalpha_nexp
    ))
  }
  
  # Compute outer model matrix
  X0 <- object$xt$basis$evalX(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros.
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0)
  
  attr(Xtot, "inner_linpred_unscaled") <- xsm_unscaled
  return(Xtot)
}
