.build_nested_bspline_basis <- function(object, data, knots, si){
  
  si$xseq <- xseq <- qnorm(1:(n-1)/n, 0, sqrt(si$vr))
  
  dsmo <- object$bs.dim
  di <- length( si$alpha )
  dtot <- dsmo + di
  
  m <- object$p.order
  if(length(m) > 2){
    m <- object$p.order <-  m[1:2]
    warning("Multiple penalties on the outer smooth are not allowed.")
  }
  
  if( is.null(knots[[object$term]]) ){
    ko <- dsmo + m[1] + 2  # Number of outer B-spline knots
    knots[[object$term]] <- 3 * qt((1:ko)/(ko+1), df = 3)
  }

  # Call this just to get info such as rank, null-space etc
  out <- .my.smooth.construct.bs.smooth.spec(object, data, knots)
  
  # Create B-spline basis
  basis <- basis_bspline(knots = knots[[object$term]], m = m[1])
  
  if(max(abs(out$X - basis$evalX(data[[object$term]])$X0)) > 1e-6){ 
    stop("Problem during B-splines basis construction") 
  }

  # Effect is not centered, so we impose that sum_i f(x_i) = 0 were x_i ~ N(0, vr).
  # We need to a) create X0 corresponding to x_i ~ N(0, vr), 
  #            b) Add a global slope
  #            c) find null space (NS) of colMeans(X0)
  #            d) project original X and S on NS
  tmpX <- basis$evalX(xseq)$X0
  tmpX <- cbind(xseq, tmpX)
  NS <- Null( colMeans(tmpX) )
  out$X <- cbind(data[[object$term]], out$X) %*% NS
  
  # Here a) bs.dim stays the same: added slope (+1) & removed intercept (-1)
  #      b) rank of pen stays the same unless full-rank (in which case must decrease by 1)
  #      c) null.space decreases by one, unless is was already empty.
  if(m[2] > 1){ 
    # Need to penalise 1st derivative too
    tmp <- object
    tmp$p.order[2] <- 1
    S1 <- .my.smooth.construct.bs.smooth.spec(tmp, data, knots)$S[[1]]
    out$S[[1]] <- out$S[[1]] + 1e-2 * S1 
    out$rank <- out$rank + 1
    out$null.space.dim <- out$null.space.dim - 1
  } else {
    # if m[2] == 1 then null space stays same: added slope & removed intercept
    if(m[2] == 0){ out$null.space.dim <- 1 } # Global slope is null space
  }
  
  # Global slope is unpenalised hence pad S with zeros and reparametrise
  out$S[[1]] <- t(NS) %*% rbind(0, cbind(0, out$S[[1]])) %*% NS 
  
  # Reparametrise the outer smooth so that penalty is diagonal
  sm <- gamFactory:::.diagPen(X = out$X, S = out$S[[1]], out$rank)
  
  # Model matrix includes inner and outer matrix (inner is just dummy) 
  out$X <- cbind(matrix(0, n, di), sm$X) 
  
  # Both penalty matrices are diagonal diag( c(0, 0, 0, ..., 1, 1, 1, ..., 0, 0)) with 
  # as many 1s as rank of penalty
  out$S <- list(rbind(cbind(matrix(0, di, di), matrix(0, di, dsmo)),
                      cbind(matrix(0, dsmo, di), sm$S)))
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
  out$xt$basis <- .wrap_nested_basis(b = basis, 
                                     P = NS %*% sm$B, 
                                     slope = TRUE)
  return(out)
}