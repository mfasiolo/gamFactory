#' 
#' Nested single-index and adaptive exponential smoothing smooth constructor
#'
#' @description
#' This function constructs a smooth term of type \code{"si_nexpsm"}, 
#' which combines an inner linear projection (single-index transformation)
#' with an outer adaptive exponential smoothing layer, followed by a 
#' B-spline basis expansion. It is designed to be used within the 
#' \pkg{mgcv} framework for flexible nested smooths.
#'
#' @details
#' The smooth is composed of three main components:
#' \itemize{
#'   \item \strong{Inner layer:} a linear projection \eqn{z_t = X_{si}^\top \alpha_{si}}, 
#'         capturing a single-index transformation of covariates.
#'   \item \strong{Middle layer:} an adaptive exponential smoothing process
#'         \eqn{s_t = \omega_t s_{t-1} + (1 - \omega_t) z_t}, 
#'         where \eqn{\omega_t = \text{sigmoid}(W_t^\top \alpha_{nexp})} 
#'         determines the dynamic smoothing weight.
#'   \item \strong{Outer layer:} a B-spline basis expansion on the smoothed
#'         output \eqn{s_t}, allowing smooth non-linear effects.
#' }
#'
#'
#' @param object A smooth specification object, typically generated 
#'               by \code{\link[mgcv]{s}} or \code{\link[mgcv]{te}} within 
#'               the \pkg{mgcv} model formula.
#' @param data A list containing the model frame components, including
#'             the design matrix and covariates used in the transformation.
#' @param knots An optional list of user-supplied knots for spline basis
#'              construction. Can be \code{NULL} (default).
#'
#' @return A smooth object containing the constructed design matrix, 
#'         penalty matrices, and additional attributes required by 
#'         the \pkg{mgcv} fitting engine.
#'
#' @importFrom Matrix rankMatrix
#' @export
#'
#'
smooth.construct.si_nexpsm.smooth.spec <- function(object, data, knots){
  si <- object$xt$si
  if(is.null(si)) si <- object$xt$si <- list()
  positive_si <- isTRUE(si$positive_si)
  
  ## 0. Split design matrix into X_si and X_nexp and times --------------------
  
  # split times
  Xall <- data[[object$term]] # Xall do not include response variable
  nms  <- colnames(Xall)
  i_t  <- which(nms == "times")            # index of 'times' column, if any
  times <- NULL                            # extract time column if present
  if (length(i_t)) {
    times <- Xall[ , i_t]
    Xall  <- Xall[ , -i_t, drop = FALSE]
  }
  si$times <- times
  
  all_columns <- ncol(Xall)
  # Xall = [X_si | X_nexp]

  if(length(i_t)){
    n_si <- si$n_si
    n_nexp <- si$n_nexp
    if (is.null(n_si)) {
      stop("`n_si` is NULL: the number of single-index covariates must be specified.", call. = FALSE)
    }
    
    if (is.null(n_nexp)) {
      stop("`n_nexp` is NULL: the number of exponential-smoothing covariates must be specified.", call. = FALSE)
    }
    
    all_columns_true <- n_si + n_nexp
    nrep <- ceiling(all_columns / all_columns_true)
    # nrep: umber of repeats
    
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
    colnames(X_si_long)   <- paste0("x", 1:n_si)      # c("East","South","West","North")
    colnames(X_nexp_long) <- paste0("w", 1:n_nexp)    # c("intercept","time_diff")
    
    # --- original X_si, X_nexp ---
    X_si   <- X_si_long
    X_nexp <- X_nexp_long
    
    # drop INF rows
    ok_si   <- rowSums(is.finite(X_si))   == ncol(X_si)
    ok_nexp <- rowSums(is.finite(X_nexp)) == ncol(X_nexp)
    
    last_good <- max(which(ok_si & ok_nexp))
    
    X_si   <- X_si[  1:last_good, , drop = FALSE]
    X_nexp <- X_nexp[1:last_good, , drop = FALSE]
    
    si$n_si <- n_si
    si$n_nexp <- n_nexp
  }
  else{
    #get n_si and n_nexp based on si,nexp attributes
    if (!is.null(si$n_nexp)){
      n_si <- all_columns - si$n_nexp
      n_nexp <- si$n_nexp
    } else{
      n_nexp <- all_columns - si$n_si
      n_si <- si$n_si
    }
    
    # check if n_si+n_nexp = all_columns
    if (n_si + n_nexp != all_columns) {
      stop(sprintf("Inconsistent number of variables: n_si + n_nexp = %d + %d != %d (ncol of Xall)", n_si, n_nexp, all_columns))
    }
    
    si$n_si <- n_si
    si$n_nexp <- n_nexp
    
    X_si <- Xall[, 1:n_si,       drop = FALSE]
    X_nexp <- Xall[, (n_si+1):(n_si+n_nexp), drop = FALSE]
  }
  
  ## 1. Optionally center X -----------------------------------------------------
  if (isTRUE(si$center)){
    X_si <- scale(X_si, scale = FALSE)
    si$xm <- attr(X_si, "scaled:center")
  } else {
    si$xm <- rep(0, n_si)
  }

  ## 2. Linear projection layer: design + penalty + initialization --------------
  S_si <- si$S_si
  noPen_si <- is.null(S_si) && is.null(si$pord)
  
  if (noPen_si){  # Case [a] no penalisation
    si$X_si    <- X_si
    si$B_si    <- diag(n_si)
    si$rank_si <- 0
  } else {
    if (is.null(S_si)){  # Case [b] "P-splines" penalty
      S_si  <- .psp(d = n_si, ord = si$pord)
      rank_si <- n_si - si$pord
    } else {  # Case [c] custom penalty Si
      rank_si <- Matrix::rankMatrix(S_si)
    }
    diag_pen_si <- gamFactory:::.diagPen(X = X_si, S = S_si, r = rank_si)
    si$S_si   <- diag_pen_si$S
    si$rank_si <- diag_pen_si$rank
    si$B_si    <- diag_pen_si$B
    si$X_si    <- diag_pen_si$X
  }
  # Initialize alpha_si and alpha_center
  if( is.null(si$alpha_center) ){
    if( noPen_si ){
      si$alpha_center <- rep(0, n_si)
    } else {
      si$alpha_center <- rep(1, n_si)
    }
  }

  if( is.null(si$alpha_si) ){ 
    if( is.null(si$alpha_center) || all(si$alpha_center == 0) ){
      si$alpha_si <- rep(1, n_si) 
    } else {
      si$alpha_si <- if (positive_si) rep(0.01, n_si) else rep(0, n_si)
    }
  }
  
  # Reparametrise and then impose that variance should be 1
  si$alpha_si  <- solve(si$B_si, si$alpha_si)
  si$alpha_center <- solve(si$B_si, si$alpha_center)
  
  tmp <- sd(si$X_si %*% (si$alpha_si + si$alpha_center))
  
  si$alpha_si      <- si$alpha_si / tmp
  si$alpha_center  <- si$alpha_center / tmp
  
  browser()
  
  if (positive_si) {
    # transform to log scale to ensure positive
    si$alpha_si <- log(pmax(si$alpha_si, 1e-8))
  }
  
  ## 3. Exponential smoothing layer: design + penalty + init α ------------------
  S_nexp <- si$S_nexp
  
  noPen_nexp <- is.null(S_nexp)
  
  if (noPen_nexp){
    si$X_nexp     <- X_nexp
    si$B_nexp     <- diag(n_nexp)
    si$rank_nexp <- 0
  } else {
    rAlp_nexp <- Matrix::rankMatrix(S_nexp)
    diag_pen_nexp <- gamFactory:::.diagPen(X = X_nexp, S = S_nexp, r = rAlp_nexp)
    si$S_nexp     <- diag_pen_nexp$S
    si$rank_nexp <- diag_pen_nexp$rank
    si$B_nexp     <- diag_pen_nexp$B
    si$X_nexp     <- diag_pen_nexp$X
  }
 
  
  # Initialize alpha_nexp without alpha_scale
  if (is.null(si$alpha_nexp)){
    si$alpha_nexp <- rep(0,n_nexp)
    g <- deriv_si_nexp(X_si = si$X_si, 
                       X_nexp = si$X_nexp, 
                       param = c(si$alpha_nexp,si$alpha_si),
                       times = times, 
                       alpha_center = si$alpha_center,
                       Z0 = si$Z0,
                       positive_si = positive_si)$d0
  } else {
    alpha_nexp <- solve(si$B_nexp) %*% si$alpha_nexp
    si$alpha_nexp <- alpha_nexp
    g <- deriv_si_nexp(X_si = si$X_si, 
                       X_nexp = si$X_nexp, 
                       param = c(si$alpha_nexp,si$alpha_si),
                       times = times, 
                       alpha_center = si$alpha_center,
                       Z0 = si$Z0,
                       positive_si = positive_si)$d0
  }
  
  ## 4. Center and scale the initialized inner linear predictor -----------------
  # alpha_si_inner has different scales
  if (positive_si) {
    si$alpha_si <- si$alpha_si - log(sd(g)) 
  } else {
    si$alpha_si <- si$alpha_si / sd(g)      
  }
  
  si$alpha_center <- si$alpha_center / sd(g) 
  si$sd_g   <- sd(g)
  
  s     <-  (g - mean(g))/sd(g)   # standardization
  
  data[[object$term]] <- s            # pass to outer smooth
  
  si$s <- s  
  si$mean_g <- mean(g)  
  si$tmp <- tmp

  ## 5. Build outer B-spline basis ----------------------------------------------
  # temporarily merge alpha_si and alpha_nexp for basis building
  si$alpha <- c(si$alpha_nexp, si$alpha_si)  
  
  out <- gamFactory:::.build_nested_bspline_basis(object = object,
                                                  data   = data,
                                                  knots  = knots,
                                                  si     = si)
  
  ## 6. Add penalties to output object ------------------------------------------
  ### 6.1 alpha_si-layer penalty (if any)
  if (!noPen_si) {
    
    di_si <- ncol(si$X_si)
    di_nexp <- ncol(si$X_nexp)
    di_beta <- out$bs.dim - di_si - di_nexp
    
    si   <- out$xt$si
    S_si <- si$S_si
    S_nexp <- matrix(0, di_nexp, di_nexp)
    S_beta <- matrix(0, di_beta, di_beta)
    
    out$S[[2]] <- as.matrix(bdiag(S_nexp, S_si, S_beta))
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank_si)
    out$rank <- c(out$rank, si$rank_si)
  }
  
  ### 6.2 α-layer penalty (α₀ unpenalized)
  if (!noPen_nexp) {
    
    di_si <- ncol(si$X_si)
    di_nexp <- ncol(si$X_nexp)
    di_beta <- out$bs.dim - di_si - di_nexp
    
    S_nexp <- si$S_nexp
    S_beta <- matrix(0, di_beta, di_beta)
    S_si <- matrix(0, di_si, di_si)
    
    out$S[[3]] <- as.matrix(bdiag(S_nexp, S_si, S_beta))
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank_nexp)
    out$rank <- c(out$rank, si$rank_nexp)
  }
  
  S2 <- if (!is.null(out$S) && length(out$S) >= 2) out$S[[2]] else NULL
  S3 <- if (!is.null(out$S) && length(out$S) >= 3) out$S[[3]] else NULL
  
  if (is.null(S2) && !is.null(S3)) {
    out$S[[2]] <- out$S[[3]]
    out$S[[3]] <- NULL
  }
  ## 7. Return ------------------------------------------------------------------
  class(out) <- c("si_nexpsm", "nested")
  return(out)
}
