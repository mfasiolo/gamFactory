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
  
  # Information on single index matrix and penalty
  si <- object$xt$si

  dsi <- ncol( si$X )
  n <- nrow( si$X )
  
  ## a truncated power spline constructor method function
  ## object$p.order = null space dimension
  if( length(object$p.order)==1 ){ object$p.order <- c(3, 2) }
  m <- object$p.order
  
  # Construct P-spline basis
  out <- smooth.construct.ps.smooth.spec(object, data, knots)
  
  dsmo <- out$bs.dim
  dtot <- dsmo + dsi
  
  # If intercept is not penalised, we penalise average value of smooth over normal distributed
  # vector of single index values 
  if( m[2] != 0 ){
    xseq <- qnorm(1:(n-1)/n, 0, object$xt$si$vr)
    tmp <- smoothCon(object = s(x, bs = "ps", k = dsmo, m = m),
                     data = data.frame(x = xseq),
                     knots = list(x = xlim), scale.penalty = FALSE)[[1]]
    xme <- colMeans( splines::spline.des(tmp$knots, x = xseq, ord = tmp$m[1] + 2, outer.ok = T)$design )
    out$rank <- out$rank + 1
    out$null.space.dim <- out$null.space.dim - 1 
    out$S[[1]] <- out$S[[1]] + xme %*% t(xme)
  }
  
  # Reparametrise the outer smooth so that penalty is diagonal
  sm <- .diagPen(X = out$X, S = out$S[[1]], out$rank)
  
  out$X <- cbind(si$X, sm$X) 
  
  # Both penalty matrices are diagonal diag( c(0, 0, 0, ..., 1, 1, 1, ..., 0, 0)) with as many 1s as rank of penalty
  if( !out$fixed ){ 
    out$S <- list(rbind(cbind(matrix(0, dsi, dsi), matrix(0, dsi, dsmo)),
                        cbind(matrix(0, dsmo, dsi), sm$S)))
    out$S[[2]] <- rbind(cbind(si$S, matrix(0, dsi, dsmo)),
                        cbind(matrix(0, dsmo, dsi), matrix(0, dsmo, dsmo)))
  }
  out$bs.dim <- dtot
  out$null.space.dim <- c(dtot - sm$rank, dtot - si$rank)
  out$rank <- c(sm$rank, si$rank)
  out$D <- NULL
  out$df <- dtot     
  out$C <- matrix(0, 0, dtot)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$updateX <- TRUE
  out$xt$splineDes <- constrSplineDes("k" = dsmo, "m" = m, "lim" = xlim, "B" = sm$B)
  
  class(out) <- "si.smooth"
  return( out )
} 
