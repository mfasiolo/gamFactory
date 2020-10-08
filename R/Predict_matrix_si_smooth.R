#'
#' Predict using single index effects
#' 
#' @name Predict.matrix.si.smooth
#' @rdname Predict.matrix.si.smooth
#' @export
#'
Predict.matrix.si.smooth <- function(object, data){
  # x-limits for P-spline basis
  xlim <- sort( object$xt$xlim )
  if( is.null(xlim) ){ xlim <- c(-6, 6) }
  
  si <- object$xt$si
  
  xa <- object$xt$si$xa
  if( is.null(xa) ){
    xa <- si$X %*% si$alpha
  } 
  
  out <- object$xt$splineDes(x = xa, deriv = 0)$X0
  out <- cbind(matrix(0, nrow(out), ncol(si$X)), out) 

  return(out)
}
