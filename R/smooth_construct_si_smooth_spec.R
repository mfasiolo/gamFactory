#'
#' Single index effects for mgcv
#' 
#' @name smooth.construct.si.smooth.spec
#' @rdname smooth.construct.si.smooth.spec
#' @importFrom MASS Null
#' @export
#'
smooth.construct.si.smooth.spec <- function(object, data, knots){
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # Inner model matrix (to be projected via single index)
  Xi <- data[[object$term]]
  
  # Need to center Xi and save colMeans because we need to subtract is when using new data
  Xi <- scale(Xi, scale = FALSE)
  si$xm <- attr(Xi, "scaled:center")
  
  di <- ncol( Xi )
  n <- nrow( Xi )
  
  # Information on single index matrix and penalty is in "si"
  # Reparametrise Xi so that the penalty on the single index vector is diagonal
  
  if( is.null(si$vr) ){ si$vr <- 1 }
  
  Si <- si$S
  if( is.null(Si) ){ 
    if( is.null(si$pord) ){ si$pord <- 1 }
    Si <- .psp(d = di, ord = si$pord)
    rankSi <- ncol(Xi) - si$pord
  } else {
    rankSi <- rankMatrix(Si)
  }
  
  si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = rankSi))
  
  # Need to initialize inner coefficient? If so, alpha chosen so that var(X %*% alpha) = si$vr 
  alpha <- si$alpha
  if( is.null(alpha) ){ alpha <- si$alpha <- rep(1, di) * sqrt(si$vr) / sd(rowSums(si$X)) }
  
  ax <- drop( si$X %*% alpha )
  data[[object$term]] <- ax
  
  out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # Add inner penalty matrix
  dsmo <- out$bs.dim - di
  si <- out$xt$si
  out$S[[2]] <- rbind(cbind(si$S, matrix(0, di, dsmo)),
                      cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo)))
  out$null.space.dim <- c(out$null.space.dim, out$bs.dim - si$rank)
  out$rank <- c(out$rank, si$rank)
  
  class(out) <- "si.smooth"
  return( out )
} 
