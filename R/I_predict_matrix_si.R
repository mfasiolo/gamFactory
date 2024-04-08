#'
#' Predict using single index effects
#' 
#' @noRd
.predict.matrix.si <- function(object, data){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  alpha <- si$alpha
  a0 <- si$a0
  
  # Need to subtract colMeans of original data "xm" and rescale using B
  Xi <- t(t(data[[object$term]]) - si$xm)  %*% si$B
  
  xa <- Xi %*% (alpha + a0)
  
  # Compute outer model matrix
  X0 <- object$xt$basis$evalX(x = xa, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xa

  return(Xtot)
}
