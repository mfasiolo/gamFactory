#'
#' Predict using nested exponential smoothing effects
#' 
#'
#' @noRd
.predict.matrix.nexpsm <- function(object, data){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  alpha <- si$alpha
  a0 <- alpha[1]
  a1 <- alpha[-1]
  
  if( is.null(si$xm) ){
    si$xm <- 0
  }
  
  # Extract variables from data frame:
  # - "y" is variable to be exponentially smoothed (stored in vector or matrix)
  # - "x" is model matrix used to model the exp smoothing rate
  # - "times" is the vector of times that at which the response variable of the 
  #           GAM is observed (hence not the same at the GAM above). 
  Xi <- data[[object$term]]
  n <- nrow( Xi )
  nms <- colnames(Xi)
  x <- as.vector( Xi[ , which(nms == "y")] )
  times <- NULL
  tmp <- which(nms == "times")
  if( length(tmp) ){
    times <- Xi[ , tmp]
  }
  Xi <- Xi[ , which(nms == "x"), drop = FALSE]
  nrep <- ceiling( length(x)/n )
  Xi <- lapply(0:(nrep-1), function(ii){
    dXi <- ncol(Xi) / nrep
    as.matrix(Xi[ , (ii*dXi + 1):(dXi*(ii+1))])
  })
  Xi <- do.call("rbind", Xi)
  
  # Need to rescale using B
  Xi <- Xi  %*% si$B
  
  xsm_unscaled <- expsmooth(y = x, Xi = Xi, beta = a1, times = times)$d0 - si$xm
  xsm <-  exp(a0) * xsm_unscaled
  
  # Compute outer model matrix
  X0 <- object$xt$basis$evalX(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xsm_unscaled
  
  return(Xtot)
  
}

