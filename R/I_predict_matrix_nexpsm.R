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
  
  # Need to rescale using B
  Xi <- data[[object$term]][ , -1, drop = FALSE]  %*% si$B
  
  xsm_unscaled <- expsmooth(y = data[[object$term]][ , 1], Xi = Xi, beta = a1)$d0 - si$xm
  xsm <-  exp(a0) * xsm_unscaled
  
  # Compute outer model matrix
  X0 <- object$xt$basis$evalX(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xsm_unscaled
  
  return(Xtot)
  
}

