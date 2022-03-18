.build_nested_pspline_basis <- function(object, data, knots, si){
  
  xlim <- sort( object$xt$xlim )  
  xseq <- si$xseq
  di <- length( si$alpha )
  
  ## a truncated power spline constructor method function
  ## object$p.order = null space dimension
  if( length(object$p.order)==1 ){ object$p.order <- c(3, 2) }
  m <- object$p.order
  
  # Construct initial P-spline basis
  out <- smooth.construct.ps.smooth.spec(object, data, knots)
  
  # Effect is not centered, so we impose that sum_i f(x_i) = 0 were x_i ~ N(0, vr).
  # We need to a) create X0 corresponding to x_i ~ N(0, vr), 
  
  #            a) Add a global slope
  #            b) find null space (NS) of colMeans(X0)
  #            c) project original X and S on NS
  tmp <- smoothCon(object = s(x, bs = "ps", k = out$bs.dim, m = m),
                   data = data.frame(x = xseq),
                   knots = list(x = xlim), scale.penalty = FALSE)[[1]]
  tmp$X <- cbind(xseq, tmp$X)
  NS <- Null( colMeans(tmp$X) )
  
  out$X <- cbind(data[[object$term]], out$X) %*% NS
  
  # Here a) bs.dim stays the same: added slope (+1) & removed incercept (-1)
  #      b) rank of pen stays the same unless full-rank (in which case must decrease by 1)
  #      c) null.space decreases by one, unless is was already empty.
  if(m[2] > 1){ 
    # Need to penalise slope in P-spline basis
    out$S[[1]] <- out$S[[1]] + 1e-2 * .Ppen(o = 1, k = out$bs.dim) 
    out$rank <- out$rank + 1
    out$null.space.dim <- out$null.space.dim - 1
  } else {
    # if m[2] == 1 then null space stays same: added slope & removed intercept
    if(m[2] == 0){ out$null.space.dim <- 1 } # Global slope is null space
  }
  
  # Global slope is unpenalised
  out$S[[1]] <- t(NS) %*% rbind(0, cbind(0, out$S[[1]])) %*% NS 
  
  dsmo <- out$bs.dim
  dtot <- dsmo + di
  
  # Reparametrise the outer smooth so that penalty is diagonal
  sm <- gamFactory:::.diagPen(X = out$X, S = out$S[[1]], out$rank)
  
  # Model matrix includes inner and outer matrix 
  out$X <- cbind(matrix(0, n, di), sm$X) 
  
  # Both penalty matrices are diagonal diag( c(0, 0, 0, ..., 1, 1, 1, ..., 0, 0)) with 
  # as many 1s as rank of penalty
  if( !out$fixed ){ 
    out$S <- list(rbind(cbind(matrix(0, di, di), matrix(0, di, dsmo)),
                        cbind(matrix(0, dsmo, di), sm$S)))
  }
  out$bs.dim <- dtot
  out$null.space.dim <- dtot - sm$rank
  out$rank <- sm$rank
  out$D <- NULL
  out$df <- dtot     
  out$C <- matrix(0, 0, dtot)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$plot.me <- FALSE
  out$repara <- TRUE # But nothing will happen as penalty is diagonal
  
  # Extra stuff needed later on. 
  out$xt$si <- si
  out$xt$basis <- .wrap_nested_basis(b = basis_bspline(k = dsmo, m = m[1], lim = xlim), 
                                     P = NS %*% sm$B, 
                                     slope = TRUE)
  return(out)
}