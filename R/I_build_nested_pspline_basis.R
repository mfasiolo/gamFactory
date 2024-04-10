# Construct outer B-spline basis for nested smooth effects
# 
.build_nested_bspline_basis <- function(object, data, knots, si){
  
  dsmo <- object$bs.dim    # Number of outer basis functions before imposing constraints (see below)
  di <- length( si$alpha ) # Number of transformation coefficients
  
  term <- object$term
  x <- data[[term]]
  knots_x <- knots[[term]] 
  
  if( abs(mean(x)) > 1e-6 ){
    stop("The data should be centered!")
  }
  
  m <- object$p.order
  if(length(m) > 2){
    m <- object$p.order <-  m[1:2]
    warning("Multiple penalties on the outer smooth are not allowed.")
  }
  
  # Determine the knots of the outer B-spline basis
  if( is.null(knots_x) ){
    ko <- dsmo + m[1] + 2      # Total number of knots B-spline knots
    kin <- ko - 2 * (m[1] + 1) # Number of inner knots
    kex <- c(-6, 6)            # Inner knots go from a to b uniformly
    dx <- diff(kex) / (kin-1)   
    knots_x <- c(seq(kex[1] - dx*(m[1]+1), kex[1]-dx, dx),  
                              seq(kex[1], kex[2], dx), 
                              seq(kex[2]+dx, kex[2] + dx*(m[1]+1), dx))
    knots[[object$term]] <- knots_x
  }
  
  # Call this just to get info such as rank, null-space etc...
  out <- .my.smooth.construct.bs.smooth.spec(object, data, knots)

  # Create B-spline basis
  basis <- basis_bspline(knots = knots_x, m = m[1])
  
  # Just checking that the model matrices correspond for observations falling
  # withing the knot range.
  tmp <- which(x >= basis$krange[1] & x <= basis$krange[2])
  if(max(abs(out$X[tmp, ] - basis$evalX(x)$X0[tmp, ])) > 1e-6){ 
    stop("Problem during B-splines basis construction") 
  }
  
  n <- nrow(out$X)
 
  # Outer smooth effect is not centered, so we impose that f(0) = 0, i.e. a point constraint at the origin.
  # We also impose, for j = 1 and 2, that f''(krange[j]) = f'''(krange[j]) = f''''(krange[j]) =  0. So in total we have 7 constraints.
  # If the order of derivative penalty >= 2, we also orthogonalise B-spline bases wrt linear basis (straight line). So 8 contraints in this case.
  # We need to a) find null space (NS) corresponding to constraints, to impose the constraints via reparametrisation
  #            b) if add_slope == TRUE we do some adjustment to "separate" linear basis from the rest 
  #            c) project original model matrix X on NS and reparametrise S
  #            d) reduce bs.dim by the number of constraints
  #            e) numerically calculate rank of S (I don't know how to get it analytically)
  X00 <- basis$evalX(x = 0)
  Xth <- basis$evalX(x = basis$krange, deriv = 4)
  con <- rbind(X00$X0, Xth$X2, Xth$X3, Xth$X4)
  NS <- Null(t(con))
  
  # We add constrain to remove linear basis from first "dsmo" columns, but then 
  # we add an unpenalised linear basis in the last column of the model matrix.
  # The NS re-parametrisation matrix is padded with 0s and (a single) 1 to make sure linear
  # effect is in last column of reparametrised X. The S penalty is padded with 0s to
  # make sure that the linear basis is unpenalised. We add to Xth the value of the
  # derivatives of the linear basis at the extreme knots.
  add_slope <- m[2] >= 2
  if(add_slope){
    con <- rbind(con, t(x)%*%out$X)
    NS <- Null(t(con))
    NS <- rbind(cbind(NS, rep(0, nrow(NS))), c(rep(0, ncol(NS)), 1))
    out$X <- cbind(out$X, x)
    out$S[[1]] <- rbind(cbind(out$S[[1]], 0), 0)
    out$bs.dim <- out$bs.dim + 1
    Xth$X0 <- cbind(Xth$X0, basis$krange)
    Xth$X1 <- cbind(Xth$X1, 1)
    Xth$X2 <- cbind(Xth$X2, 0)
    Xth$X3 <- cbind(Xth$X3, 0)
    Xth$X4 <- cbind(Xth$X4, 0)
  }
  out$X <- out$X %*% NS
  out$S[[1]] <- t(NS) %*% out$S[[1]] %*% NS 
  ncon <- nrow(NS) - ncol(NS)
  out$bs.dim <- out$bs.dim - ncon 
  dsmo <- out$bs.dim
  out$rank <- rankMatrix(out$S[[1]])
  
  # Linearly extrapolate outer B-spline basis beyond inner knots (those at which we imposed the constraints above)
  # NOTE method = "smooth" won't work if add_slope == TRUE
  out$X <- linextr(x = x, b = list(X0 = out$X), th = basis$krange, 
                   Xbo = Xth$X0%*%NS, Xbo1 = Xth$X1%*%NS, method = "simple")$X0 # If you change "simple" here, you need to do it also elsewhere!!
  
  # To aid identifiability we make penalty on the outer smooth full rank
  # if(out$rank < dsmo){
  #   eig <- eigen(out$S[[1]])
  #   eig$values[-(1:out$rank)] <- pmin(eig$values[1] * 1e-5, eig$values[out$rank])
  #   out$S[[1]] <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  #   out$rank <- ncol(out$X)
  # }
  
  # Reparametrise the outer smooth so that penalty is diagonal (we need diagonal penalty otherwise mgcv will
  # reparametrise again with Sl.repara)
  sm <- gamFactory:::.diagPen(X = out$X, S = out$S[[1]], out$rank)
  
  # Final model matrix includes inner and outer matrix (inner part first "di" columns, just dummy needed by mgcv) 
  out$X <- cbind(matrix(0, n, di), sm$X) 
  
  # Pad penalty on outer smooth with zeros corresponding to the parameters of the inner transformation.
  # Penalty on inner coefficients will be added outside this function. 
  # Both inner and outer penalty matrices are diagonal diag( c(0, 0, 0, ..., 1, 1, 1, ..., 0, 0)) with 
  # as many 1s as the rank of penalty.
  out$S <- list(rbind(cbind(matrix(0, di, di), matrix(0, di, dsmo)),
                      cbind(matrix(0, dsmo, di), sm$S)))
  dtot <- dsmo + di
  out$bs.dim <- dtot
  out$null.space.dim <- dtot - sm$rank
  out$rank <- sm$rank
  out$D <- NULL
  out$df <- dtot     
  out$C <- matrix(0, 0, dtot)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$plot.me <- FALSE
  out$repara <- TRUE # But mgcv will not reparametrised at penalty is diagonal 
                     # (but, I think, setting this to FALSE leads to an error)
  
  # Store extra things needed later on. 
  out$xt$si <- si
  out$xt$basis <- .wrap_nested_basis(b = basis, P = NS %*% sm$B, Xth = Xth, add_slope = add_slope)
  return(out)
}