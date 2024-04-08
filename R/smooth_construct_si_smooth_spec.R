#'
#' Single index effects for mgcv
#' 
#' @name smooth.construct.si.smooth.spec
#' @rdname smooth.construct.si.smooth.spec
#' @importFrom MASS Null
#' @export
#'
smooth.construct.si.smooth.spec <- function(object, data, knots){
  
  # Most information on single index matrix and penalty is inside "si" list.
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # Inner model matrix (to be projected via single index)
  Xi <- data[[object$term]]
  
  # Need to center Xi and save colMeans because we need to subtract it when predicting using new data
  Xi <- scale(Xi, scale = FALSE)
  si$xm <- attr(Xi, "scaled:center")
  
  di <- ncol( Xi )
  n <- nrow( Xi )
  
  # Dealing with inner penalty
  Si <- si$S
  no_pen <- is.null(Si) && is.null(si$pord)
  if( no_pen ){ # Case [a] no penalisation
    si$X <- Xi
    si$B <- diag(nrow = ncol(Xi))
    si$rank <- 0 
  } else {
    if( is.null(Si) ){ # Case [b] "P-splines" penalty
      if( is.null(si$pord) ){ si$pord <- 0 }
      Si <- .psp(d = di, ord = si$pord)
      rankSi <- ncol(Xi) - si$pord
    } else { # Case [c] custom penalty Si
      rankSi <- rankMatrix(Si)
    }
    # Reparametrise Xi so that the penalty on the single index vector is diagonal
    si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = rankSi))
  }
  
  # alpha is vector of inner coefficients, si$alpha is a vector of initial values for it.
  # alpha0 is an offset such that the full_alpha = alpha + alpha0
  if( is.null(si$a0) ){
    if( no_pen ){
      si$a0 <- rep(0, di)
    } else {
      si$a0 <- rep(1, di)
    }
  }
  if( is.null(si$alpha) ){ 
    if( is.null(si$a0) || all(si$a0 == 0) ){
      si$alpha <- rep(1, di) 
    } else {
      si$alpha <- rep(0, di) 
    }
  }
  # Reparametrise and then impose that variance should be 1
  si$alpha <- solve(si$B, si$alpha)
  si$a0 <- solve(si$B, si$a0)
  tmp <- sd(si$X %*% (si$alpha + si$a0))
  si$alpha <- si$alpha / tmp
  si$a0 <- si$a0 / tmp
  
  # Compute single index vector and store it in the data
  ax <- drop( si$X %*% (si$alpha + si$a0) ) 
  data[[object$term]] <- ax
  
  # Construct the B-splines corresponding to the outer smooth effect 
  out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # Add inner penalty matrix (diagonalised and padded with zeros corresponding to the outer coefficients)
  if( !no_pen ){
    dsmo <- out$bs.dim - di
    si <- out$xt$si
    out$S[[2]] <- rbind(cbind(si$S, matrix(0, di, dsmo)),
                        cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo)))
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank)
    out$rank <- c(out$rank, si$rank)
  }
  
  class(out) <- c("si", "nested")
  return( out )
} 
