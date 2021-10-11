#'
#' Nested MGKS effects for mgcv
#' 
#' @name smooth.construct.mgks.smooth.spec
#' @rdname smooth.construct.mgks.smooth.spec
#' @examples 
#' n <- 100
#' n0 <- 50
#' p <- 5
#' dat <- data.frame(y = rnorm(n))
#' dat$Xi <- matrix(rnorm(n*p), n, p)
#' X0 <- matrix(rnorm(n0*p), n0, p)
#' x <- rnorm(n0)
#' aaa <- smoothCon(s(Xi, bs = "mgks", m = c(3, 2),
#'             xt = list(si = list(X0 = X0, x = x), sumConv = FALSE)),
#'           data = dat)
#' @export
smooth.construct.mgks.smooth.spec <- function(object, data, knots)
{ 
  si <- object$xt$si
  
  Xi <- data[[object$term]]
  di <- (ncol(Xi)-1)/2 + 1
  n <- nrow( Xi )
  n0 <- si$n0
  si$x <- x <- Xi[1:n0, 1]
  si$X0 <- X0 <- Xi[1:n0, 2:di, drop = FALSE]
  si$X <- Xi <- Xi[ , -(1:di), drop = FALSE]
  
  if( is.null(si$vr) ){ si$vr <- var(x) }
  
  # x-limits for outer P-spline basis
  xlim <- sort( object$xt$xlim )
  if( is.null(xlim) ){ xlim <- range(x) }
  
  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = si$vr (target variance)
    tmp <- mgks(y = x, X = Xi, X0 = X0, beta = -log(colSds(X0)))$d0
    alpha <- si$alpha <- c(log(sqrt(si$vr)/sd(tmp)), -log(colSds(X0))) 
    data[[object$term]] <- exp(alpha[1]) * tmp
  } else {
    data[[object$term]] <- exp(alpha[1]) * mgks(y = x, X = Xi, X0 = X0, beta = alpha[-1])$d0
  }

  ## A truncated power spline constructor method function
  ## object$p.order = null space dimension
  if( length(object$p.order)==1 ){ object$p.order <- c(3, 2) }
  m <- object$p.order
  
  # Construct initial P-spline basis
  out <- smooth.construct.ps.smooth.spec(object, data, knots)
  
  # Effect is not centered, so we impose that sum_i f(x_i) = 0 were x_i is raw data (no smoothing).
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
  #      c) null.space decreases by one, unless it was already empty.
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
  }
  out$bs.dim <- dtot
  out$null.space.dim <- c(dtot - sm$rank)
  out$rank <- c(sm$rank)
  out$D <- NULL
  out$df <- dtot     
  out$C <- matrix(0, 0, dtot)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$plot.me <- FALSE
  out$repara <- TRUE # But nothing will happen as penalty is diagonal
  
  # Extra stuff needed later on. 
  # NB: "k" = dsmo+1 because we lost 1 dimension via centering constraint
  out$xt$si <- si
  out$xt$basis <- basis_pspline(k = dsmo+1, m = m[1], lim = xlim, P = NS %*% sm$B)

  class(out) <- "mgks.smooth"
  return( out )
} 
