#'
#' Predict using single index effects
#' 
#' @name Predict.matrix.si.smooth
#' @rdname Predict.matrix.si.smooth
#' @export Predict.matrix.si.smooth
#'
Predict.matrix.si.smooth <- function(object, data){
  # x-limits for P-spline basis
  xlim <- sort( object$xt$xlim )
  if( is.null(xlim) ){ xlim <- c(-6, 6) }
  
  X <- object$xt$si$X

  out <- Predict.matrix.pspline.smooth(object, data)
  out <- cbind(X[1:nrow(out), ], out)
  
  return(out)
}