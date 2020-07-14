#'
#' Single index effects for mgcv
#' 
#' @name smooth.construct.si.smooth.spec
#' @rdname smooth.construct.si.smooth.spec
#' @export smooth.construct.si.smooth.spec
#'
smooth.construct.si.smooth.spec <- function(object, data, knots)
{ 
  # Extra args:
  # - object$xt$X (n x p) matrix of variables to be projected via single index model
  # - object$xt$xlim limit knots for P-spline basis
  # Secret storage in output object:
  # - object$X stores single index matrix in the last p columns. First p columns are useless.
  # - object$xlim store limit know position
  # - object$S is (k+p) x (k+p) and the last p x p submatrix on the bottom-right is zero.
  #   Hence no penalty on the single index coefficients.
  # x-limits for P-spline basis
  xlim <- sort( object$xt$xlim )
  if( is.null(xlim) ){ xlim <- c(-6, 6) }
  
  XI <- object$xt$X
  dsi <- ncol(XI)
  n <- nrow(XI)
  
  ## a truncated power spline constructor method function
  ## object$p.order = null space dimension
  if( length(object$p.order)==1 ){ object$p.order <- c(3, 2) }
  m <- object$p.order
  
  out <- smooth.construct.ps.smooth.spec(object, data, knots)
  
  dsmo <- out$bs.dim
  
  B <- solve( t( sqrtPenDec(P = out$S[[1]], r = out$rank) ))
  
  # print(out$S[[1]])
  # print(D)
  
  # Reparametrise so that P-spline penalty is diagonal
  out$X <- out$X %*% B 
  out$X <- cbind(XI, out$X) 
  
  # Expand S to include 0 penalty on SI coefficients
  if( !out$fixed ){ 
    out$S <- list(rbind(cbind(matrix(0, dsi, dsi), matrix(0, dsi, dsmo)),
                        cbind(matrix(0, dsmo, dsi), diag(1, dsmo))))
  }
  out$bs.dim <- out$bs.dim + dsi
  out$null.space.dim <- dsi
  out$rank <- dsmo
  
  out$df <- ncol(out$X)     
  out$C <- matrix(0, 0, ncol(out$X))
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$D <- NULL
  
  class(out)<-"si.smooth"  # Give object a class
  return( out )
} 
