#'
#' Predict using single index effects
#' 
#' @name Predict.matrix.si.smooth
#' @rdname Predict.matrix.si.smooth
#' @export
#'
Predict.matrix.si.smooth <- function(object, data){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  alpha <- si$alpha
  
  d <- length( alpha )
  
  # Need to subtract colMeans of original data "xm" and rescale using B
  Xi <- t(t(matrix(data[[object$term]], ncol = d)) - si$xm)  %*% si$B
  
  xa <- Xi %*% alpha 
  
  # Compute outer model matrix
  X0 <- object$xt$splineDes(x = xa, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), d), X0) 

  return(Xtot)
}
