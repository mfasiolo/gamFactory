#'
#' Single index effects for mgcv
#' 
#' @name smooth.construct.si.smooth.spec
#' @rdname smooth.construct.si.smooth.spec
#' @importFrom MASS Null
#' @export
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

  # Inner model matrix (to be projected via single index)
  Xi <- data[[object$term]]
  
  # Need to center Xi and save colMeans because we need to subtract is when using new data
  Xi <- scale(Xi, scale = FALSE)
  object$xt$si$xm <- attr(Xi, "scaled:center")
  
  dsi <- ncol( Xi )
  n <- nrow( Xi )
  
  # Information on single index matrix and penalty is in "si"
  # Reparametrise Xi so that the penalty on the single index vector is diagonal
  si <- object$xt$si
  if( is.null(si$vr) ){ si$vr <- 1 }
  if( is.null(si$ord) ){ si$ord <- 1 }
  si <- append(si, gamFactory:::.diagPen(X = Xi, S = .psp(d = dsi, ord = si$ord), r = ncol(Xi) - si$ord))
  
  # Need to initialize inner coefficient? If so, alpha chosen so that var(X %*% alpha) = si$vr 
  alpha <- si$alpha
  if( is.null(alpha) ){ alpha <- si$alpha <- rep(1, dsi) * sqrt(si$vr) / sd(rowSums(si$X)) }
  
  ax <- si$X %*% alpha
  data[[object$term]] <- ax
  
  ## a truncated power spline constructor method function
  ## object$p.order = null space dimension
  if( length(object$p.order)==1 ){ object$p.order <- c(3, 2) }
  m <- object$p.order
  
  # Construct initial P-spline basis
  out <- smooth.construct.ps.smooth.spec(object, data, knots)
  
  # Effect is not centered, so we impose that sum_i f(x_i) = 0 were x_i ~ N(0, vr).
  # We need to a) create X0 corresponding to x_i ~ N(0, vr), 
  #            b) find null space (NS) of xme = colMeans(X0)
  #            c) project original X and S on NS
  xseq <- qnorm(1:(n-1)/n, 0, sqrt(si$vr))
  tmp <- smoothCon(object = s(x, bs = "ps", k = out$bs.dim, m = m),
                   data = data.frame(x = xseq),
                   knots = list(x = xlim), scale.penalty = FALSE)[[1]]
  xme <- colMeans( splines::spline.des(tmp$knots, x = xseq, ord = tmp$m[1] + 2, outer.ok = T)$design )
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
  dtot <- dsmo + dsi
  
  # Reparametrise the outer smooth so that penalty is diagonal
  sm <- gamFactory:::.diagPen(X = out$X, S = out$S[[1]], out$rank)
  
  # Model matrix includes inner and outer matrix 
  out$X <- cbind(matrix(0, n, dsi), sm$X) 
  
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
  out$plot.me <- FALSE
  
  # Extra stuff needed later on. 
  # NB: "k" = dsmo+1 because we lost 1 dimension via centering constraint
  out$xt$si <- si
  out$xt$splineDes <- PsplineDesign(k = dsmo+1, m = m[1], lim = xlim, P = NS %*% sm$B)
  
  class(out) <- "si.smooth"
  return( out )
} 
