.my.smooth.construct.bs.smooth.spec <- function(object,data,knots) {
  # a B-spline constructor method function
  # Here m[1] is as in Wood's 2016 book section 5.3.3 (i.e., m = 2 for cubic)
  # but we need to increase it by 1 to make Predict.matrix.Bspline.smooth happy later on....
  m <- object$p.order  # m[1] - basis order, m[2] - penalty order
  m[1] <- m[1] + 1     
  object$m <- object$p.order <- m # NOTE object$m and $p.order modified, will fix at the end
  
  x <- data[[object$term]]  
  k <- knots[[object$term]]
  
  # Total number of knots, number of interior knots
  nk <- length(k)
  nkin <- nk - 2 * m[1]

  object$X <- splines::spline.des(k, x, m[1]+1, x*0, outer.ok = TRUE)$design # get model matrix
  
  if (length(unique(x)) < object$bs.dim) warning("basis dimension is larger than number of unique covariates")
  
  ## now construct derivative based penalty. Order of derivate
  ## is equal to m, which is only a conventional spline in the 
  ## cubic case...
  object$knots <- k; 
  class(object) <- "Bspline.smooth"  # Give object a class
  k0 <- k[m[1] + 1:nkin] ## the interior knots
  object$D <- object$S <- list()
  m2 <- m[2:length(m)] ## penalty orders
  if (length(unique(m2))<length(m2)) stop("multiple penalties of the same order is silly")
  for (i in 1:length(m2)) { ## loop through penalties
    object$deriv <- m2[i] ## derivative order of current penalty
    pord <- m[1]-m2[i] ## order of derivative polynomial 0 is step function
    if (pord<0) stop("requested non-existent derivative in B-spline penalty") 
    h <- diff(k0) ## the difference sequence...
    ## now create the sequence at which to obtain derivatives
    if (pord==0) k1 <- (k0[2:nkin]+k0[1:(nkin-1)])/2 else {
      h1 <- rep(h/pord,each=pord)
      k1 <- cumsum(c(k0[1],h1)) 
    } 
    dat <- data.frame(k1);names(dat) <- object$term 
    D <- Predict.matrix.Bspline.smooth(object,dat) ## evaluate basis for mth derivative at the k1
    object$deriv <- NULL ## reset or the smooth object will be set to evaluate derivs in prediction! 
    if (pord==0) { ## integrand is just a step function...
      object$D[[i]] <- sqrt(h)*D
    } else { ## integrand is a piecewise polynomial...
      P <- solve(matrix(rep(seq(-1,1,length=pord+1),pord+1)^rep(0:pord,each=pord+1),pord+1,pord+1))
      i1 <- rep(1:(pord+1),pord+1)+rep(1:(pord+1),each=pord+1) ## i + j
      H <- matrix((1+(-1)^(i1-2))/(i1-1),pord+1,pord+1)
      W1 <- t(P)%*%H%*%P
      h <- h/2 ## because we map integration interval to to [-1,1] for maximum stability
      ## Create the non-zero diagonals of the W matrix... 
      ld0 <- rep(sdiag(W1),length(h))*rep(h,each=pord+1)
      i1 <- c(rep(1:pord,length(h)) + rep(0:(length(h)-1) * (pord+1),each=pord),length(ld0))
      ld <- ld0[i1] ## extract elements for leading diagonal
      i0 <- 1:(length(h)-1)*pord+1
      i2 <- 1:(length(h)-1)*(pord+1)
      ld[i0] <- ld[i0] + ld0[i2] ## add on extra parts for overlap
      B <- matrix(0,pord+1,length(ld))
      B[1,] <- ld
      for (k in 1:pord) { ## create the other diagonals...
        diwk <- sdiag(W1,k) ## kth diagonal of W1
        ind <- 1:(length(ld)-k)
        B[k+1,ind] <- (rep(h,each=pord)*rep(c(diwk,rep(0,k-1)),length(h)))[ind]  
      }
      ## ... now B contains the non-zero diagonals of W
      B <- bandchol(B) ## the banded cholesky factor.
      ## Pre-Multiply D by the Cholesky factor...
      D1 <- B[1,]*D
      for (k in 1:pord) {
        ind <- 1:(nrow(D)-k)
        D1[ind,] <- D1[ind,] + B[k+1,ind] * D[ind+k,]
      }
      object$D[[i]] <- D1
    }
    object$S[[i]] <- crossprod(object$D[[i]])
  }
  object$rank <- object$bs.dim-m2  # penalty rank 
  object$null.space.dim <- min(m2)    # dimension of unpenalized space 
  object$interior.knots.bound <- k0[c(1, length(k0))]
  object$m[1] <- object$m[1] - 1 ###### Need to reset these, See above
  object$p.order[1] <- object$p.order[1] - 1 

  return( object )
} ### end of B-spline constructor