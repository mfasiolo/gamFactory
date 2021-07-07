#'
#' Predict using nested exponential smoothing effects
#' 
#' @name Predict.matrix.nexpsm.smooth
#' @rdname Predict.matrix.nexpsm.smooth
#' @export
#'
Predict.matrix.nexpsm.smooth <- function(object, data){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  alpha <- si$alpha
  
  # Need to rescale using B
  Xi <- data[[object$term]][ , -1, drop = FALSE]  %*% si$B
  
  xsm <-  expSmooth(y = data[[object$term]][ , 1], Xi = Xi, beta = alpha)$d0 
  
  # Compute outer model matrix
  X0 <- object$xt$splineDes(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  return(Xtot)
  
}
