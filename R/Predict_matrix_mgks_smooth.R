#'
#' Predict using MGKS smooth effects
#' 
#' @name Predict.matrix.mgks.smooth
#' @rdname Predict.matrix.mgks.smooth
#' @export
#'
Predict.matrix.mgks.smooth <- function(object, data){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  alpha <- si$alpha
  
  Xi <- data[[object$term]]
  di <- (ncol(Xi)-1)/2
  n0 <- si$n0

  xsm <- mgks(y = Xi[1:n0, 1], X = Xi[ , -(1:(di+1)), drop = FALSE], 
              X0 = Xi[1:n0, 2:(di+1), drop = FALSE], beta = alpha)$d0 
  
  # Compute outer model matrix
  X1 <- object$xt$splineDes(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X1), di), X1) 
  
  return(Xtot)
  
}
