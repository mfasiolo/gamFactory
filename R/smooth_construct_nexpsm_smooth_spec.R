#'
#' Nested adaptive exponential smoothing effects for mgcv
#' 
#' @name smooth.construct.nexpsm.smooth.spec
#' @rdname smooth.construct.nexpsm.smooth.spec
#' @export
#'
smooth.construct.nexpsm.smooth.spec <- function(object, data, knots)
{ 
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # Covariate to be exponentially smoothed
  # Inner model matrix and penalty (used to model exponential smoothing coefficient) 
  Xi <- data[[object$term]]
  x <- Xi[ , 1]
  Xi <- Xi[ , -1]
  n <- length( x )
  Si <- si$Si
  di <- ncol(Xi) + 1
  
  if( is.null(si$vr) ){ si$vr <- var(x) }
  
  # x-limits for outer P-spline basis
  if( is.null(object$xt$xlim) ){ object$xt$xlim <- c(-6, 6) * sqrt(si$vr) }
  
  # Reparametrise Xi so that the penalty on alpha is diagonal
  si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = si$rankp))

  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = si$vr (target variance)
    g <- expsmooth(y = x, Xi = si$X, beta = rep(0, di-1))$d0
    alpha <- si$alpha <- c(log(sqrt(si$vr)/sd(tmp)), rep(0, di-1))
  } else {
    g <- expsmooth(y = x, Xi = si$X, beta = alpha[-1])$d0
  }
  
  # Center and scale the initialiized inner linear preditor
  data[[object$term]] <- exp(alpha[1]) * (g - mean(g))
  
  si$x <- x
  si$xseq <- qnorm(1:(n-1)/n, 0, sqrt(si$vr))
  
  out <- .build_nested_pspline_basis(object = object, data = data, knots = knots, si = si)
  
  # Add inner penalty matrix
  dsmo <- out$bs.dim - di
  si <- out$xt$si
  out$S[[2]] <- rbind(cbind(rbind(0, cbind(0, si$S)), matrix(0, di, dsmo)),
                      cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo)))
  out$null.space.dim <- c(out$null.space.dim, out$bs.dim - si$rank)
  out$rank <- c(out$rank, si$rank)
  
  class(out) <- "nexpsm.smooth"
  return( out )
} 
