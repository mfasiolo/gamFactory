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
  xlim <- sort( object$xt$xlim )
  if( is.null(xlim) ){ xlim <- range(x) }
  
  # Reparametrise Xi so that the penalty on alpha
  si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = si$rankp))

  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = si$vr (target variance)
    tmp <- expsmooth(y = x, Xi = si$X, beta = rep(0, di-1))$d0
    alpha <- si$alpha <- c(log(sqrt(si$vr)/sd(tmp)), rep(0, di-1)) 
    data[[object$term]] <- exp(alpha[1]) * tmp
  } else {
    data[[object$term]] <- exp(alpha[1]) * expsmooth(y = x, Xi = si$X, beta = alpha[-1])$d0
  }
  
  ## A truncated power spline constructor method function
  ## object$p.order = null space dimension
  if( length(object$p.order)==1 ){ object$p.order <- c(3, 2) }
  m <- object$p.order
  
  # Construct initial P-spline basis
  out <- smooth.construct.ps.smooth.spec(object, data, knots)
  
  # Effect is not centered, so we impose that sum_i f(x_i) = 0 were x_i is raw data (no exp smoothing).
  # We need to a) find null space (NS) of xme = colMeans(X0)
  #            b) project original X and S on NS
  tmp <- smoothCon(object = s(x, bs = "ps", k = out$bs.dim, m = m),
                   data = data.frame(x = x),
                   knots = list(x = xlim), scale.penalty = FALSE)[[1]]
  xme <- colMeans( splines::spline.des(tmp$knots, x = x, ord = tmp$m[1] + 2, outer.ok = T)$design )
  NS <- Null( xme %*% t(xme) )
  out$X <- out$X %*% NS
  out$S[[1]] <- t(NS) %*% out$S[[1]] %*% NS 
  
  # Here a) bs.dim certainly decreases by 1 
  #      b) rank of pen stays the same unless full-rank (in which case must decrease by 1)
  #      c) null.space decreases by one, unless is was already empty.
  out$bs.dim <- out$bs.dim - 1
  out$rank <- min(out$rank, out$bs.dim)
  out$null.space.dim <- max(out$null.space.dim - 1, 0) 
  
  dsmo <- out$bs.dim
  dtot <- dsmo + di
  
  # Reparametrise the outer smooth so that penalty is diagonal
  sm <- gamFactory:::.diagPen(X = out$X, S = out$S[[1]], out$rank)
  
  # Model matrix includes inner and outer matrix 
  out$X <- cbind(matrix(0, n, di), sm$X) 
  
  # Both penalty matrices are diagonal diag( c(0, 0, 0, ..., 1, 1, 1, ..., 0, 0)) with as many 1s as rank of penalty
  if( !out$fixed ){ 
    out$S <- list(rbind(cbind(matrix(0, di, di), matrix(0, di, dsmo)),
                        cbind(matrix(0, dsmo, di), sm$S)))
    out$S[[2]] <- rbind(cbind(rbind(0, cbind(0, si$S)), matrix(0, di, dsmo)),
                        cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo)))
  }
  out$bs.dim <- dtot
  out$null.space.dim <- c(dtot - sm$rank, dtot - si$rank)
  out$rank <- c(sm$rank, si$rank)
  out$D <- NULL
  out$df <- dtot     
  out$C <- matrix(0, 0, dtot)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$plot.me <- FALSE
  out$repara <- FALSE
  
  # Extra stuff needed later on. 
  # NB: "k" = dsmo+1 because we lost 1 dimension via centering constraint
  out$xt$si <- si
  out$xt$basis <- basis_pspline(k = dsmo+1, m = m[1], lim = xlim, P = NS %*% sm$B)
  out$xt$si$x <- x

  class(out) <- "nexpsm.smooth"
  return( out )
} 
