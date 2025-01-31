#'
#' Predict using nested exponential smoothing effects
#' 
#'
#' @noRd
.predict.matrix.nexpsm <- function(object, data, get.xa = FALSE){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  # Do NOT remove this! We need it in postproc_gam_nl.
  if( is.null(si$xm) ){
    si$xm <- 0
  }
  
  alpha <- si$alpha
  a0 <- alpha[1]
  a1 <- alpha[-1]
  
  # Extract variables from data frame:
  # - "y" is variable to be exponentially smoothed (stored in vector or matrix)
  # - "x" is model matrix used to model the exp smoothing rate
  # - "times" is the vector of times that at which the response variable of the 
  #           GAM is observed (hence not the same at the GAM above). 
  Xi <- data[[object$term]]
  n <- nrow( Xi )
  nms <- colnames(Xi)
  x <- as.vector( t(Xi[ , which(nms == "y")]) )
  times <- NULL
  tmp <- which(nms == "times")
  if( length(tmp) ){
    times <- Xi[ , tmp]
  }
  Xi <- Xi[ , which(nms == "x"), drop = FALSE]
  nrep <- ceiling( length(x)/n )
  dXi <- ncol(Xi)/nrep
  if(nrep > 1){
    tmp <- rep(1:dXi, nrep)
    Xi <- apply(Xi, 1, function(x) do.call("cbind", tapply(x, tmp, I)), simplify = FALSE)
    Xi <- do.call("rbind", Xi)
  }
  
  # Need to rescale using B
  Xi <- Xi  %*% si$B
  
  xsm_list <- expsmooth(y = x, Xi = Xi, beta = a1, times = times, deriv = get.xa)
  xsm_unscaled <- xsm_list$d0 - si$xm
  xsm <-  exp(a0) * xsm_unscaled
  
  if(get.xa){ 
    return(list(xa = xsm,
                xa_da = exp(a0)*xsm_list$d1))
  }
  
  # Compute outer model matrix
  X0 <- object$xt$basis$evalX(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xsm_unscaled
  
  return(Xtot)
  
}

