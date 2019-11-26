## (c) Simon N Wood 2019 Provided under GPL2

gam.fit1 <- function(G,Sl,sp,scale,Sld=NULL,start=NULL,etastart=NULL,mustart=NULL,control = gam.control(),nt=1,stab=TRUE,ident=FALSE) {
## Cholesky based gam fit for gradient only based estimation via REML.
## G is a pre-fit gam object. Sl is a penalty object created by Sl.setup(G).
## * sp contains smoothing parameters. Last element may be log scale parameter if this is not fixed.
##   for extended family the first G$family$n.theta elements are the extra parameters, theta.
## * Before first calling this:
##   1. any extended family pre-initialization must be called
##   2. Sl must be created using Sl.setup. This latter involves a reparameterization.
##   3. The re-parameterization from 2 must be applied to the model matrix G$X.
##   4. For exponential families the family should have been fixed.
##   5. Assumption is that L is handled externally.
##   6. For discrete methods the Sl reparameterization has to be undone before
##      using discrete matrix multiplication routines, as these have to operate on raw parameterizattion
## note: currently no gamma parameter
## note: should indef XWX be handled via fisher weights - currently pivoted chol is simply truncated...?
## NOTE: currently assumes start is re-parameterized: but could this clash with any setting of start in pre-initialize code??
##       -- probably only an issue for discrete methods, as compact model matrix can not be reparameterized in advance.

  wtoc <- function(gamma) {
    ## transform from working parameters to regression coefficients
    if (is.null(gamma)) return(gamma)
    if (!is.null(G$g.index) && any(G$g.index)) gamma[G$g.index] <- exp(gamma[G$g.index])
    gamma
  } ## wtoc

  ctow <- function(beta) {
    ## transform from regression coefficients to working paramters
    if (!is.null(G$g.index) && any(G$g.index)) beta[G$g.index] <- log(beta[G$g.index])
    beta
  } ## ctow

  raw <- function(beta) {
    ## get raw regression coefs from working...
    if (discrete) { ## must be raw parameterization if discrete
      raw.beta <- mgcv:::Sl.repa(ddetS$rp,beta,l=-1)
      raw.beta <- wtoc(raw.beta)
      if (length(idrop)) { ## insert 0 for unidentifiable
        rb <- raw.beta; raw.beta <- rep(0,length(rb)+length(idrop))
	raw.beta[-idrop] <- rb
      }
      raw.beta <- Sl.initial.repara(Sl,raw.beta,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1)
    } else raw.beta <- wtoc(beta)
    raw.beta
  } ## raw

  if (is.null(Sld)) {
    Sld <- Sl
    idrop <- rep(0,0)
  } else {
    idrop <- attr(Sld,"drop") ## index of terms to drop in discrete case, because unidentifiable 
  }

  y <- G$y
  discrete <-  !is.null(G$Xd) ## is the model matrix discretized?
  nop <- ncol(G$X) ## number of parameters before any that still need to be dropped for strucural reasons (discrete case only) 
  x <- if (discrete) list(Xd = G$Xd,lpid=G$lpid,kd = G$kd,ks=G$ks,ts=G$ts,dt=G$dt,qc=G$qc,v=G$v,p=nop) else G$X
  if (discrete) nop <- nop - length(idrop) ## number of identifiable parameters
  if (discrete) {
    attr(x,"lpi") <- attr(G$X,"lpi")
    if (inherits(G$family,"general.family")&&is.null(G$family$discrete.ok)) stop("family does not support discrete computation methods")
  }
  mu <- NULL
  nobs <- length(y)
  weights <- G$w
  family <- G$family
  n.theta <- 0

  if (inherits(G$family,"extended.family")) if (inherits(G$family,"general.family")) {
    efam <- FALSE;gfam <- TRUE
    llf <- G$family$ll
    E <- attr(Sl,"E") ## E'E is a regularizing penalty which may be used in initialization
    ## for discrete methods E needs to be back in original parameterization
    if (discrete) E <- Sl.initial.repara(Sl,E,inverse=TRUE,both.sides=FALSE,cov=FALSE) else
    if (length(idrop)) E <- E[,-idrop,drop=FALSE] ## make sure this is consistent with X
    attr(E,"use.unscaled") <- TRUE ## signal initialization code that E not to be further scaled
  } else { ## preinitialize extended family
    efam <- TRUE;gfam <- FALSE
    if (family$n.theta>0) { ## there are extra parameters to estimate
      n.theta <- family$n.theta
      ind <- 1:family$n.theta
      theta <- sp[ind] ## parameters of the family
      family$putTheta(theta)
      sp <- sp[-ind]   ## log smoothing parameters
    } else theta <- family$getTheta() ## fixed value
  } else gfam <- efam <- FALSE
  
  if (scale>0) scale.known <- TRUE else {
    ## unknown scale parameter, trial value supplied as 
    ## final element of sp. 
    scale.known <- FALSE
    nsp <- length(sp)
    scale <- exp(sp[nsp])
    sp <- sp[-nsp]
  }
  rho <- sp

  offset <- G$offset
  
  linkinv <- family$linkinv
  dev.resids <- family$dev.resids
  if (!efam && !gfam) {
      variance <- family$variance
      mu.eta <- family$mu.eta
      if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
  }
  
  ## run initialization code for family...
  valideta <- family$valideta
  if (is.null(valideta))
        valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu))
        validmu <- function(mu) TRUE
	
  start0 <- start;start <- NULL	
  if (is.null(mustart)) {
        eval(family$initialize)
  } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
  }
  ## initialization is in raw parameterization for discrete... need to convert...
  if (discrete && !is.null(start)) {
    start <- Sl.initial.repara(Sl,start,inverse=FALSE,both.sides=FALSE,cov=FALSE,nt=1)
    if (length(idrop)) start <- start[-idrop]
  }

  ## if initialization produces 'start' and some parameters are +ve only
  if (!is.null(G$g.index) && !is.null(start)) start[G$g.index] <- log(pmax(1e-6,start[G$g.index])) 

  if (!is.null(start0)) start <- start0 ## prevent initialization from over-writing supplied 'start'
  
  if (!gfam&&is.matrix(y)&&ncol(y)>1) stop("This family should not have a matrix response")

  ## get derivatives of log|S|_+ and update Sl with smoothing parameters...
  fixed <- rep(FALSE,length(rho))

  ddetS <- mgcv:::ldetS(Sld,rho,fixed,np=nop,root=TRUE,repara=stab,nt=nt,deriv=1)
  
  Sld <- ddetS$Sl ## has smoothing parameters updated and reparameterized penalties in Srp
  St <- crossprod(ddetS$E) ## total penalty matrix (reparameterized)
    
  if (is.null(start)) { ## never true for general families
    eta <- if (!is.null(etastart)) drop(etastart) else family$linkfun(drop(mustart))
    bSb <- 0
  } else {
    #raw.beta <- if (discrete) Sl.initial.repara(Sl,start,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1) else start ## start must be in raw param for discrete
    if (discrete && length(idrop)) { ## need to re-insert zeroes for non-identified for discrete 
      raw.beta <- rep(0,length(start)+length(idrop))
      raw.beta[-idrop] <- wtoc(start)
    } else raw.beta <- wtoc(start)
    if (discrete) raw.beta <- Sl.initial.repara(Sl,raw.beta,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1) ## back to raw parameterization
    ## repara neither (or both) of x & start before this... 
    if (!gfam) eta <- offset + if (discrete) mgcv:::Xbd(G$Xd,raw.beta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else drop(x %*% raw.beta) 
    start <- mgcv:::Sl.repara(ddetS$rp,start,inverse=FALSE,both.sides=FALSE) ## apply stabilizing reparameterization
    Sb <- mgcv:::Sl.mult(Sld,start)
    bSb <- sum(start*Sb)
  }

  if (!discrete) x <- mgcv:::Sl.repara(ddetS$rp,x,inverse=FALSE,both.sides=FALSE) ## apply stabilizing reparameterization

  if (!discrete) raw.beta <- wtoc(start)
  beta <- start

  if (gfam) {
    dcheck <- FALSE
    ll <- llf(y,x,raw.beta,weights,G$family,offset=G$offset,deriv=as.numeric(dcheck))
    if (dcheck) { ## derivative testing of llf
      eps <- 1e-4
      fd.lb <- ll$lb
      fd.lbb <- ll$lbb
      for (i in 1:length(raw.beta)) {
        rb1 <- raw.beta;rb1[i] <- rb1[i] + eps
	ll1 <- llf(y,x,rb1,weights,G$family,offset=G$offset,deriv=1)
	fd.lb[i] <- (ll1$l-ll$l)/eps
	fd.lbb[,i] <- (ll1$lb-ll$lb)/eps
      }
      plot(fd.lb,ll$lb)
      plot(fd.lbb,ll$lbb)
    }
    ddev <- dev <- -2 * ll$l + bSb
  } else {
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
         stop("cannot find valid starting values: please specify some")
    ddev <- dev <- if (is.null(start)) 1.1*sum(dev.resids(y, family$linkinv(0) , weights)) else sum(dev.resids(y, mean(mu), weights)) + bSb
    ## ... avoid converging or wild divergence at iter 1
  }
  conv <- FALSE
  strictly.additive <- family$family=="gaussian" && family$link=="identity" && is.null(G$g.index)

  fisher <- FALSE ## note: unused at present, but what if XWX n +ve def?

  if (efam) theta <- family$getTheta()

  ## The current penalty matrix is required
  steepest <- FALSE

  ## in what follows, raw.beta is regression coefs, in raw parameterization for discrete,
  ## beta is working coefficients, may have been re-parameterized and involve non-linear
  ## parameterization. penalties apply to beta, model matrix to raw.beta.

  if (!is.null(G$g.index)&&is.null(start)) {
    ## need first fit iteration in regression parameterization
    ## to get starting values.
    g.index <- G$g.index;G$g.index <- NULL 
  } else g.index <- NULL

  eps <- 1e-4
  ####################################################
  for (iter in 1:control$maxit) { ## main fitting loop

    dcheck <- FALSE
    K <- if (dcheck) length(beta) + 1 else 1
    for (k in 1:K) { ## FD loop
    if (k>1) {
      beta <- beta0;beta[k-1] <- beta[k-1] + eps
      raw.beta <- raw(beta)
      if (!gfam) {
        eta <- offset + if (discrete) mgcv:::Xbd(G$Xd,raw.beta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else drop(x %*% raw.beta)
        mu <- linkinv(eta)
      }	
    }
    if (gfam) {
      ll <- llf(y,x,raw.beta,weights,G$family,offset=G$offset,deriv=1)
      g <- -ll$lb
      XWXS <- -ll$lbb
    } else if (efam) {
      dd <- mgcv:::dDeta(y,mu,weights,theta,family,0) ## derivatives of deviance w.r.t. eta
      w <- dd$Deta2 * .5  ## so X'diag(w)X is Hessian of deviance
      g <- if (discrete) mgcv:::XWyd(G$Xd,rep(.5,nobs),dd$Deta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else
                                              drop(dd$Deta %*% x)*.5 ## grad of deviance wrt coefs
      XWXS <- if (discrete) mgcv:::XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nt,G$drop) else crossprod(x,w*x) ## unpenalized Hessian
    } else { ## exponential family
      
      var.val <- variance(mu)
      if (any(is.na(var.val))) stop("NAs in V(mu)")  
      if (any(var.val == 0)) stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")
            
      good <- (weights > 0) & (mu.eta.val != 0)
         
      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("No observations informative at iteration %d", iter))
        break
      }

      if (fisher) { ## Conventional Fisher scoring
	w <- weights * mu.eta.val^2/var.val 
      } else { ## full Newton
      	alpha <- 1+ (y-mu)*(family$dvar(mu)/var.val + family$d2link(mu)*mu.eta.val)
        alpha[alpha==0] <- .Machine$double.eps
        w <- alpha * weights * mu.eta.val^2/var.val 
      }
     
      
      g <- if (discrete) -mgcv:::XWyd(G$Xd,mu.eta.val/var.val,y-mu,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else
                         -drop(((y-mu)*mu.eta.val/var.val) %*% x)  ## grad of deviance wrt coefs
      XWXS <- if (discrete) mgcv:::XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nt,G$drop) else crossprod(x,w*x) ## unpenalized Hessian
    }
   
    ## deal with any non-linearity of parameters...
    if (!is.null(G$g.index)) {
      XWX0 <- XWXS ## store for log det deriv calc
      g0 <- g
      bg <- exp(beta[G$g.index])
      g[G$g.index] <- g[G$g.index]*bg
      XWXS[G$g.index,] <- bg* XWXS[G$g.index,]
      XWXS[,G$g.index] <- t(bg* t(XWXS[,G$g.index]))
      XWXS[G$g.index,G$g.index] <- XWXS[G$g.index,G$g.index] + diag(g[G$g.index]) 
    } else XWX0 <- NULL ## so grad and hess are now w.r.t. beta
    XWX <- XWXS ## save here to return at end
    
    if (discrete) { ## have to impose reparameterization (as it can't be pre-absorbed)
      XWXS <- Sl.initial.repara(Sl,XWXS,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
      g <- Sl.initial.repara(Sl,g,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
      if (length(idrop)) {
        XWXS <- XWXS[-idrop,-idrop]
	g <- g[-idrop]
      }	
      XWXS <- mgcv:::Sl.repa(ddetS$rp,XWXS,l=-2,r=-1) 
      g <- mgcv:::Sl.repa(ddetS$rp,g,r=-1)
    }
    ## Now add penalty to the Hessian...
    if (!steepest) {
      XWXS <- XWXS + St
      ## form pivoted Cholesky decomposition with diagonal pre-conditioning...
      d <- diag(XWXS);ind <- d<=0
      d[ind] <- 1;d[!ind] <- sqrt(d[!ind])
      R <- if (nt>1) mgcv:::pchol(t(XWXS/d)/d,nt) else suppressWarnings(chol(t(XWXS/d)/d,pivot=TRUE))
      r <- min(attr(R,"rank"),Rrank(R))
      p <- ncol(XWXS)
      piv <- attr(R,"pivot") #;rp[rp] <- 1:p
      if (r<p) { ## drop rank deficient terms...
        R <- R[1:r,1:r]
        dpiv <- piv[(r+1):p] ## index of dropped terms 
        piv <- piv[1:r]
      }
 
      if (iter == 1 && is.null(start)) { 
	## Idea is that the initial step would be to (X'WX + S)^{-1} X'W z, where z = h'(mu)(y-mu) + eta (h is link here)
	## but the grad g = X'(y-mu)/(h'(mu)V(mu)) while W^{-1} = V(mu)h'mu^2, so -X'Wz = g - X'Weta
	## Hence if we set initial beta to zero and g to as below then initial step will be correct.
	## Also holds for extended. Not for general, but then start is never null.
        xwe <- if (discrete) -mgcv:::XWyd(G$Xd,w,eta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else -drop((w*eta) %*% x)
        if (discrete) {
	  xwe <- Sl.initial.repara(Sl,xwe,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
	  if (length(idrop)) xwe <- xwe[-idrop]
	  xwe <- mgcv:::Sl.repa(ddetS$rp,xwe,l=-2) 
        }
        beta <- rep(0,p)
	raw.beta <- rep(0,p+length(idrop))
        g <- g + xwe
      }
      ## now complete gradient by adding penalty deriv
    } ## !steepest
    Sb <- mgcv:::Sl.mult(Sld,beta)
    g <- g + Sb
    if (dcheck) {
      if (k==1) {
        beta0 <- beta
        gb <- g
	XWXf <- XWXS*0
      } else {
        XWXf[,k-1] <- (g-gb)/eps
      }

    }
    } ## FD loop end

    step <- rep(0,p)
    devtol <- control$epsilon*abs(dev)
    ## test convergence (crucial to do this before either branch - was on Newton branch only, leading to
    ## a perfect zig-zag on exact rank deficiency)

    
    
    
    
    
    
    
    ###### CHANGE FROM
    # conv <- (ddev <= dev * control$epsilon && !any(abs(g)>devtol))
    ###### TO
    # print(iter)
    # print("max abs gradient:")
    # print(max(abs(g)))
    # print("ddev:")
    # print(ddev)
    # print("devtol:")
    # print(devtol)
    conv <- (ddev <= devtol && !any(abs(g)>devtol))
    
    
    
    
    
    
    
    
    
    if (conv) break
    if (steepest) { ## take a steepest descent step, only in space dropped at last Newton update
      dstep <- -g[dpiv]
      dstep <- .1*(mean(abs(beta[dpiv]))+mean(abs(step))+mean(abs(beta))*1e-3)*dstep/mean(abs(dstep))
      step[dpiv] <- dstep ## equivalent to adding a +ve def perturbation to inverse Hessian, so still a descent direction
    } else { 
      ## compute the Newton step...
      step[piv] <- -backsolve(R,(forwardsolve(t(R),g[piv]/d[piv])))/d[piv]
    }
    ## try the step, halving if needed...
    beta1 <- beta + step
#    if (discrete) { ## must be raw parameterization if discrete
#      raw.beta1 <- mgcv:::Sl.repa(ddetS$rp,beta1,l=-1)
#      raw.beta1 <- wtoc(raw.beta1)
#      if (length(idrop)) { ## insert 0 for unidentifiable
#        rb <- raw.beta1; raw.beta1 <- rep(0,length(rb)+length(idrop))
#	raw.beta1[-idrop] <- rb
#      }
#      raw.beta1 <- Sl.initial.repara(Sl,raw.beta1,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1)
#    } else raw.beta1 <- wtoc(beta1)
    raw.beta1 <- raw(beta1)
    if (!gfam) {
      if (discrete) {
        eta1 <- offset + mgcv:::Xbd(G$Xd,raw.beta1,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop)
      } else eta1 <- offset +  drop(x %*% raw.beta1)
      mu1 <- linkinv(eta1)
    }  
    Sb1 <- mgcv:::Sl.mult(Sld,beta1)
   
    for (k in 1:20) { ## step control loop...
      dev.up <- dev1 <- if (gfam) -2*llf(y,x,raw.beta1,weights,G$family,offset=G$offset,deriv=0)$l else sum(dev.resids(y, mu1, weights))
      dev1 <- dev1 + sum(beta1*Sb1)
      if (is.finite(dev1) && dev1-dev <= devtol) { ## fine
        beta <- beta1;Sb <- Sb1;
	raw.beta <- raw.beta1 
	if (!gfam) { eta <- eta1; mu <- mu1}
	break 
      } else { ## need to step half
        Sb1 <- (Sb1 + Sb)/2
        beta1 <- (beta1 + beta)/2
	if (!is.null(G$g.index)) { ## raw.beta non-linear in beta
	  raw.beta1 <- raw(beta1)
#          if (discrete) {
#           raw.beta1 <- mgcv:::Sl.repa(ddetS$rp,beta1,l=-1) ## need to undo stabilizing reparameterization
#	    raw.beta1[G$g.index] <- exp(raw.beta1[G$g.index]) ## do non-linear transform
#	    if (length(idrop)) { ## insert 0 for unidentifiable ## deal with dropped
#              rb <- raw.beta1; raw.beta1 <- rep(0,length(rb)+length(idrop))
#	      raw.beta1[-idrop] <- rb
#            }
	    ## and finally undo the initial reparameterization
#	    raw.beta1 <- Sl.initial.repara(Sl,raw.beta1,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1)
#          } else { ## not discrete - so only the non-linearity to deal with
#            raw.beta1 <- beta1
#	    raw.beta1[G$g.index] <- exp(raw.beta1[G$g.index]) ## do non-linear transform
#          }
        } else raw.beta1 <- (raw.beta1 + raw.beta)/2 ## raw.beta linear in beta 
	if (!gfam) {
	  if (!is.null(G$g.index)) { ## non-linear in beta
              eta1 <- offset + if (discrete) mgcv:::Xbd(G$Xd,raw.beta1,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else drop(x %*% raw.beta1)
          } else eta1 <- (eta1 + eta)/2 ## linear in beta 
	  mu1 <- linkinv(eta1)
	}  
      }
    } ## step halving loop
    if (iter==1 && !is.null(g.index)) {
      ## first iteration operated on regression coefficient scale
      ## to get starting values. Need to get working scale starting values
      ## and corresponding dev...
      beta[g.index] <- log(pmax(1e-6,beta[g.index]))
      G$g.index <- g.index;g.index <- NULL ## restore G$g.index to signal non-linearity
      #if (discrete) {
      #      raw.beta <- mgcv:::Sl.repa(ddetS$rp,beta,l=-1) ## need to undo stabilizing reparameterization
#	    raw.beta[G$g.index] <- exp(raw.beta[G$g.index]) ## do non-linear transform
#	    if (length(idrop)) { ## insert 0 for unidentifiable ## deal with dropped
 #             rb <- raw.beta; raw.beta <- rep(0,length(rb)+length(idrop))
#	      raw.beta[-idrop] <- rb
 #           }
#	    ## and finally undo the initial reparameterization
#	    raw.beta <- Sl.initial.repara(Sl,raw.beta,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1)
 #     } else { ## not discrete - so only the non-linearity to deal with
  #          raw.beta <- beta
#	    raw.beta[G$g.index] <- exp(raw.beta[G$g.index]) ## do non-linear transform
 #     }
      raw.beta <- raw(beta)
      if (!gfam) { 
          eta <- offset + if (discrete) mgcv:::Xbd(G$Xd,raw.beta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else drop(x %*% raw.beta)
	  mu <- linkinv(eta)
      }  
      ## get the deviance corrsponding to beta in range...
      dev1 <- if (gfam) -2*llf(y,x,raw.beta,weights,G$family,offset=G$offset,deriv=0)$l else sum(dev.resids(y, mu, weights))
      Sb1 <- mgcv:::Sl.mult(Sld,beta)
      dev1 <- dev1 + sum(Sb1*beta)
    }
    if (!is.finite(dev1)||dev1-dev > devtol) stop("step halving failure")
    ddev <- abs(dev - dev1)
    dev <- dev1
    if (steepest) steepest <- FALSE else {
      if (r<p) steepest <- TRUE
    }
    if (strictly.additive) {
      conv <- TRUE; break
    }
  } ## main fitting loop
  
  if (ident) {
    if (stab==TRUE) warning("identifiability being checked with stabilizing transform on!")
    if (!conv) warning("PIRLS not converged at identifiability check stage!")
    if (r<p) return(dpiv) else return(rep(0,0))
  }
  
  if (!conv) warning("PIRLS not converged")

  if (!is.null(G$g.index) && discrete) { ## have to impose reparameterization on XWX0 and g0 before deriv calc
      XWX00 <- XWX0 ## retain this as one correction term needs this parameterization
      XWX0 <- Sl.initial.repara(Sl,XWX0,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
      g0 <- Sl.initial.repara(Sl,g0,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
      if (length(idrop)) {
        XWX0 <- XWX0[-idrop,-idrop]
	g0 <- g0[-idrop]
      }	
      XWX0 <- mgcv:::Sl.repa(ddetS$rp,XWX0,l=-2,r=-1) 
      g0 <- mgcv:::Sl.repa(ddetS$rp,g0,r=-1)
  } else XWX00 <- NULL
  
  ## get explicit inverse of penalized Hessian
  PP <- matrix(0,p,p)
  if (nt>1) {
    P <- mgcv:::pbsi(R,nt=nt,copy=TRUE) ## invert R 
    PP[piv,piv] <-  mgcv:::pRRt(P,nt) ## PP'
  } else PP[piv,piv] <- chol2inv(R)
  PP.raw <- PP <- t(PP/d)/d ## explicit inverse of penalized Hessian (unpivoted)
  
  ## Now perform implicit differentiation, use Sl.iftChol + extra for theta
  id <- mgcv:::Sl.iftChol(Sld,XX=NULL,R=R,d=d,beta=beta,piv=piv,nt=nt)
  if (efam&&family$n.theta>0) { ## implicit differentation for the theta parameters
    ## key quantity is second partial of deviance w.r.t. \beta and rho
    dd <- mgcv:::dDeta(y,mu,weights,theta,family,1) ## derivatives of deviance w.r.t. eta (to 3rd order)
    if (discrete) {
      dd.dbt <- mgcv:::XWyd(G$Xd,rep(1,nobs),dd$Detath,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop)
      ##dd.dbt <- Sl.initial.repara(Sl,dd.dbt,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
      dd.dbt <- Sl.inirep(Sl,dd.dbt,l=2)
      if (length(idrop)) dd.dbt <- if (family$n.theta>1) dd.dbt[-idrop,,drop=FALSE] else dd.dbt[-idrop]
      dd.dbt <- mgcv:::Sl.repa(ddetS$rp,dd.dbt,l=-2)
      #dd.dbt <- Sl.inirep(Sl,dd.dbt,l=1,r=0,nt=nt)
    } else dd.dbt <- t(t(dd$Detath) %*% x) ## note worth avoiding t(X)
    db <- cbind(-PP%*%dd.dbt/2,id$db)
    id$bSb1 <- c(rep(0,family$n.theta),id$bSb1) ## bSb has no direct dependence on theta
    ddetS$ldet1 <- c(rep(0,family$n.theta),ddetS$ldet1)
  } else db <- id$db

  db.raw <- db
  ## beta = exp(gamma) => dbeta/drho = dgamma/drho * dbeta/dgamma =  dgamma/drho * beta
  if (!is.null(G$g.index)) db.raw[G$g.index,] <- db.raw[G$g.index,] * raw.beta[G$g.index]

  if (discrete) { ## reversing reparameterizations
    db.raw <- mgcv:::Sl.repa(ddetS$rp,db.raw,l=-1)
    if (length(idrop)) {
      rb <- db.raw; db.raw <- matrix(0,nrow(db.raw)+length(idrop),ncol(db.raw))
      db.raw[-idrop,] <- rb
    }
    db.raw <- Sl.inirep(Sl,db.raw,l=1,r=0,nt=nt)
  } 
  
  ## compute LAML first order derivatives
  ## -\lambda_j b'S_jb/2 + dlog|S|_+/drho_j/2 - tr(Hp^{-1}(dH/d\rho + \lambda_j S_j))/2
  ## first is in id$bSb1, second in ddetS$ldet1 

  ## now the derivatives of log|X'WX+S|, tr(Hp^{-1}S_j) terms...
  
  dH <- mgcv:::d.detXXS(Sld,PP,nt=1,deriv=1) ## note: could be made more efficient given only 1st derivs

  ## Now the tr(Hp^{-1}dH/drho_j) = a'dw/drho_j terms...
  if (gfam) { ## general families
    ## call llf with deriv = 2 and fh set to inverse penalized hessian to get tr(Hp^{-1}dH/drho_j) vector
    #if (discrete) stop("discrete general not coded")
    if (!is.null(G$g.index)) { ## need to pass some extra information to llf for computing log det derivs
      if (length(idrop)) { ## pad out g.index for passing to llf
        gind <- rep(FALSE,length(idrop)+length(G$g.index))
	gind[!idrop] <- G$g.index
	attr(db.raw,"g.index") <- gind
      } else attr(db.raw,"g.index") <- G$g.index
      attr(db.raw,"beta") <- raw.beta
    } ## end of extra llf info
    
    if (discrete) { ## reverse repara
      PP.raw <- mgcv:::Sl.repa(ddetS$rp,PP,l=-1,r=-2)
      if (length(idrop)) { ## insert 0 rows and cols
        rb <- PP.raw;PP.raw <- matrix(0,length(idrop)+nrow(PP.raw),ncol(PP.raw)+length(idrop))
	PP.raw[-idrop,-idrop] <- rb
      }
      PP.raw <- Sl.initial.repara(Sl,PP.raw,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=nt) 
    } else PP.raw <- PP
    ## NOTE: problem here with discrete. raw.beta and db.raw code shared with non-gfam discrete => appears correct
    ##       PP.raw or gamlss.gH seem most probable as source of error...
    ddetH <- dH$d1 - llf(y,x,raw.beta,weights,G$family,offset=G$offset,d1b=db.raw,deriv=2,fh=PP.raw)$d1H
    lsat <- 0
  } else { ## extended and exponential families
    ## work through computation of dW/drho_j...
    if (discrete) {
      #db.raw <- Sl.inirep(Sl,db,l=1,r=0,nt=nt)
      deta <- mgcv:::Xbd(G$Xd,db.raw,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) 
    } else deta <- x %*% db.raw 
 
    if (efam) { ## extended family branch
      dW <- .5 * dd$Deta3 * deta
      if (family$n.theta) {
        dW[,1:family$n.theta] <- dW[,1:family$n.theta] + dd$Deta2th*.5
	## update part of partial of log|H| depending directly on s.p.s  
        dH$d1 <- c(rep(0,family$n.theta),dH$d1) ## no dependence on theta via penalties in this term
      }
    } else { ## exponential family
      mu.eta.val <- mu.eta(eta)
      dmu <- deta * mu.eta.val
      g1 <- 1/mu.eta.val
      g2 <- family$d2link(mu)
      g3 <- family$d3link(mu)

      V <- family$variance(mu)
      V1 <- family$dvar(mu)
      V2 <- family$d2var(mu)
      alpha1 <- -(V1/V+g2/g1) + (y-mu)*(V2/V-(V1/V)^2+g3/g1-(g2/g1)^2)
      dW <- ((alpha1 - alpha*(2*g2/g1+V1/V))/(g1^2*V))*dmu
    } ## exp fam branch

    dhc <- rep(0,ncol(dW)) ## correction to detH derivatives 
    if (!is.null(G$g.index)) {
      BHi <- PP
      bg <- exp(beta[G$g.index]) ## dbeta/dgamma
      BHi[G$g.index,] <- bg*BHi[G$g.index,]
      BHiB <- BHi
      BHiB[,G$g.index] <- t(t(BHi[,G$g.index])*bg)
      if (!discrete) a <- rowSums((x%*%BHiB)*x)
      ## now the corrections relating to the dependence of the
      ## transformation on the parameters... 
      for (k in 1:length(dhc)) {
        ## trace correction terms
	bk <- bg*db[G$g.index,k]
        BHiBk <- t(bk*t(BHi[,G$g.index]))
	## diagonal correction terms...
	XWX2 <- XWX0[,G$g.index]
	dk <- if (is.null(XWX00)) XWX0[G$g.index,]%*%db.raw[,k] else XWX00[G$g.index,]%*%db.raw[,k];
	dk <- dk*bg #dk[G$g.index] <- dk[G$g.index]*bg
	dk1 <- g0[G$g.index]; dk1 <- dk1 * bk  # dk1[G$g.index] <- dk1[G$g.index]*bk
	dk <- dk + dk1
	dhc[k] <- 2*sum(BHiBk*XWX2) + sum(diag(PP)[G$g.index]*dk)
      }
    } else {
      if (!discrete) a <- rowSums((x%*%PP)*x)
      BHiB <- PP
    }
    dW <- dW * weights ## incorporate prior weights
    dW[!is.finite(dW)] <- 0 ## drop any non-finite.
    if (discrete) { ## more efficient to explicitly work through sps 
      ddetH <- rep(0,ncol(dW))
      for (i in 1:length(ddetH)) { ## loop over smoothing/hyper parameters
         dHi <- mgcv:::XWXd(G$Xd,dW[,i],G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nt,G$drop)
	 ## apply repara...
	 dHi <- Sl.initial.repara(Sl,dHi,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
	 if (length(idrop)) dHi <- dHi[-idrop,-idrop]
	 dHi <- mgcv:::Sl.repa(ddetS$rp,dHi,l=-2,r=-1)
	 ddetH[i] <- sum(BHiB*dHi) + dH$d1[i]
      }
    } else ddetH <- colSums(a*dW) + dH$d1 ## tr(Hp^{-1}(dH/d\rho + \lambda_j S_j)
    ddetH <- ddetH + dhc
    ls <- if (efam) family$ls(y,weights,theta,scale) else family$ls(y,weights,n,scale)
    lsat <- if (efam) ls$ls else ls[1]
  } 
  laml1 <- -0.5*id$bSb1/scale + 0.5*ddetS$ldet1 - 0.5*ddetH
  if (efam&&family$n.theta>0) { ## add the direct dependence on theta to derivatives
    ind <- 1:family$n.theta
    laml1[ind] <- laml1[ind] + ls$lsth1[ind] - colSums(as.matrix(dd$Dth))/(2*scale)
  }
  if (!scale.known && !gfam) { ## append derivatives w.r.t. log scale
    ## important terms are l - b'Sb/(2*scale) - M/2 log(scale)
    ## where M is unpenalized space rank and l = l_sat - deviance/(2*scale)
    ## recall in the following that dev is penalized deviance, not raw deviance.
    dscale <- -G$min.edf/2 + 0.5*dev/scale
    ## note in next line that some idiot has coded the exponential familes to
    ## return derivs w.r.t. scale and extended w.r.t. log(scale)...
    dscale <- dscale + if (efam) ls$lsth1[length(ls$lsth1)] else ls[2]*scale
    laml1 <- c(laml1,dscale)
  } ## finished derivatives wrt log scale

  ## At least for derivative checking purposes, the LAML should also be returned.
  detH <- 2*sum(log(diag(R)*d[piv]))
  ## note that the stabilizing transform cancels between the  log|S| and the log|hess| terms
  laml <- lsat - dev/(2*scale) + ddetS$ldetS/2 - detH/2 -sum(G$min.edf)*(log(scale)-log(2*pi))/2
  dVkk <- t(db) %*% XWXS %*% db ## curvature checking matrix 

  if (inherits(G$family,"extended.family")) pearson <- NULL else {
    pearson <- sum(weights*(y-mu)^2/family$variance(mu))
    s.bar = max(-.9,mean(family$dvar(mu)*(y-mu)*sqrt(weights)/family$variance(mu)))
    if (is.finite(s.bar)) pearson <- pearson/(1+s.bar) ## divide by residual edf to get Fletcher scale estimate
  }

  aic.model <- if (inherits(G$family,"extended.family")) {
    if (inherits(G$family,"general.family")) -2*ll$l else G$family$aic(y, mu, theta, weights, dev.up)
  } else G$family$aic(y, n, mu, weights, dev.up)

  if (inherits(G$family,"general.family")) { ## get fitted values and lienar predictors.
    lpi <- attr(x,"lpi")
    if (is.null(lpi)) { ## only one...
      eta <- offset + if (discrete) mgcv:::Xbd(G$Xd,raw.beta,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else drop(x%*%beta)
      mu <- family$linkinv(eta) 
    } else { ## multiple...
      mu <- eta <- matrix(0,nobs,length(lpi))
      if (!is.null(offset)) offset[[length(lpi)+1]] <- 0
      for (j in 1:length(lpi)) {
        if (discrete) { ## could be made more efficient using lt argument of Xbd
	  beta0 <- rep(0,length(raw.beta));beta0[lpi[[j]]] <- raw.beta[lpi[[j]]]
	}  
        eta[,j] <- if (discrete) mgcv:::Xbd(G$Xd,beta0,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else as.numeric(x[,lpi[[j]],drop=FALSE] %*% beta[lpi[[j]]])
        if (!is.null(offset[[j]])) eta[,j] <-  eta[,j] + offset[[j]]
        mu[,j] <- family$linfo[[j]]$linkinv( eta[,j]) 
      }
    }
    dev.up <- NULL
  } 

  ## following 2 reparameterizations are needed for derivative checking, but
  ## could be applied in gam.fit1.postproc otherwise...
  beta <- mgcv:::Sl.repa(ddetS$rp,beta,l=-1) ## undo stabilizing reparameterization
  db <- mgcv:::Sl.repa(ddetS$rp,db,l=-1) ## undo stabilizing reparameterization

  list(laml=-laml,laml1=-laml1,fitted.values=mu,coefficients=beta,scale.est=scale,pearson=pearson,
       beta=beta,db=db,dVkk = dVkk, detH = detH,ddetH = ddetH,mu=mu,XWXS=XWXS,XWX=XWX,XWX0=XWX0,aic=aic.model,
       rp=ddetS$rp,Vb = PP.raw,discrete=discrete,iter=iter,gfam=gfam,scale.known=scale.known,n.theta=n.theta,
       n=nobs,fitted.values = mu,linear.predictors=eta,deviance=dev.up,g.index=G$g.index)
       ## NOTE: debugging returns need cleaning out eventually
} ## gam.fit1

gam.fit1.postproc <- function(b,Sl,S,off,V.rho=NULL,L=NULL,lsp0=NULL,lsp=NULL,nt=1) {
## post process a gam.fit1 fit to undo all remaining reparameterizations and produce
## Vc/b/e, edf1/2/3, F the DoF matrix and R the unpenalized XWX factor.
## b is the object returned from gam.fit1 while Sl is the penalty object sent to gam.fit1.
## S and off are smooth list and offset from initial setup.
## V.rho is approximate cov matrix for log smoothng parameters.
## NOTE: need to deal with re-insertion of any dropped parameters here, I think
  #coef <- mgcv:::Sl.repa(b$rp,b$beta,l=-1) ## undo stabilizing reparameterization
  rank <- length(b$beta)
  idrop <- attr(Sl,"drop")
  if (!is.null(b$g.index)) b$beta[b$g.index] <- exp(b$beta[b$g.index])
  if (length(idrop)) {
    np <- length(b$beta) + length(idrop)
    coef <- rep(0,np);coef[-idrop] <- b$beta
    b$beta <- coef
    if (!is.null(b$g.index)) {
      g.index <- rep(FALSE,np)
      g.index[-idrop] <- b$g.index
    }  
  } else {
    g.index <- b$g.index
    np <- length(b$beta)
  }  
  coef <- Sl.initial.repara(Sl,b$beta,inverse=TRUE,both.sides=FALSE,cov=FALSE,nt=1)
  if (b$discrete) {
    XWX <- b$XWX ## save in raw parameterization
    ## transform to fit parameterization
    b$XWX <- Sl.initial.repara(Sl,b$XWX,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt)
    if (length(idrop)) b$XWX <- b$XWX[-idrop,-idrop]
    b$XWX <- mgcv:::Sl.repa(b$rp,b$XWX,l=-2,r=-1)
  } else { ## get XWX in raw parameterization...
    XWX <- mgcv:::Sl.repa(b$rp,b$XWX,l=2,r=1)
    if (length(idrop)) {
      V <- matrix(0,np,np)
      V[-idrop,-idrop] <- XWX
      XWX <- V
    }
    XWX <-  Sl.initial.repara(Sl,XWX,inverse=TRUE,both.sides=TRUE,cov=FALSE)
  }
  
  if (b$discrete && b$gfam) Vb <- b$Vb else {
      Vb <- mgcv:::Sl.repa(b$rp,b$Vb,l=-1,r=-2)
      if (length(idrop)) {
        V <- matrix(0,np,np)
	V[-idrop,-idrop] <- Vb
	Vb <- V
      }
      Vb <- Sl.initial.repara(Sl,Vb,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=nt) 
  }
  ## need to get R, the factor of XWX
  if (!is.null(b$XWX0)) b$XWX <- b$XWX0 ## this relates to regression, not working, coefficients
  d <- sqrt(diag(b$XWX) + norm(b$XWX)*.Machine$double.eps^.5)
  R <- if (nt>1) mgcv:::pchol(t(b$XWX/d)/d,nt) else suppressWarnings(chol(t(b$XWX/d)/d,pivot=TRUE))
  ipiv <- piv <- attr(R,"pivot")
  ipiv[piv] <- 1:length(piv)
  R <- t(t(R[,ipiv])*d)
  R <- mgcv:::Sl.repa(b$rp,R,r=1) ## NOTE: gam.fit5.post.proc looks wrong in comparison
  if (length(idrop)) {
    V <- matrix(0,np,np)
    V[-idrop,-idrop] <- R
    R <- V
  }
  R <-  Sl.initial.repara(Sl,R,inverse=TRUE,both.sides=FALSE,cov=FALSE)
  ## XWX <- crossprod(R)

  ## NOTE: temporary testing inefficiency here 
  b$XWX <- mgcv:::Sl.repa(b$rp,b$XWX,l=2,r=1)
  if (length(idrop)) {
    V <- matrix(0,np,np)
    V[-idrop,-idrop] <- b$XWX
    b$XWX <- V
  }
  b$XWX <-  Sl.initial.repara(Sl,b$XWX,inverse=TRUE,both.sides=TRUE,cov=FALSE)
  ## END testing

  F <- Vb%*%XWX ## DoF matrix
  Vb <- Vb * b$scale.est
  Ve <- F%*%Vb ## 'frequentist' cov matrix
  if (!is.null(V.rho)) { ## sp uncertainty correction
    if (length(idrop)) {
      V <- matrix(0,np,ncol(b$db))
      V[-idrop,] <- b$db
      b$db <- V
    }
    db <- Sl.inirep(Sl,b$db,l=1)
    if (!is.null(L)) {
      db <- if (b$scale.known) db %*% L else db %*%  L[-nrow(L),-ncol(L)]  ## transform to derivs w.r.t. working
    }
    Vc <- Vb + if (b$scale.known) db %*% V.rho %*% t(db) else db %*% V.rho[-nrow(V.rho),-nrow(V.rho)] %*% t(db)
    Vc <- Vc + mgcv:::Vb.corr(R,L,lsp0=lsp0,S=S,off=off,dw=NULL,w=NULL,rho=lsp,Vr=V.rho,nth=b$n.theta,scale.est=!b$scale.known)*b$scale.est
  }
  edf <- diag(F)
  edf1 <- 2*edf - rowSums(t(F)*F)
  edf2 <- rowSums(Vc*XWX)
  if (!is.null(b$g.index)) { ## delta method transform all to new space...
    Vc[g.index,] <- coef[g.index]*Vc[g.index,]
    Vc[,g.index] <- t(coef[g.index]*t(Vc[,g.index]))
    Vb[g.index,] <- coef[g.index]*Vb[g.index,]
    Vb[,g.index] <- t(coef[g.index]*t(Vb[,g.index]))
    Ve[g.index,] <- coef[g.index]*Ve[g.index,]
    Ve[,g.index] <- t(coef[g.index]*t(Ve[,g.index]))
  }
  scale.est <- if (is.null(b$pearson)) b$scale.est else b$pearson/(b$n - sum(edf))
  ## NOTE: R not in correct parameterization under non-linear parameterization
  list(coef=coef,Vb=Vb,Vc=Vc,Ve=Ve,edf=edf,edf1=edf1,edf2=edf2,R=R,F=F,scale.est =scale.est,rank=rank)
} ## gam.fit1.post.proc



go.bfgs <-  function(lsp,G,Sl,scale,lsp0=NULL,Sld=NULL,L=NULL,start=NULL,etastart=NULL,
                     mustart=NULL,control = gam.control(),conv.tol=1e-6,maxNstep=2,
		     stab=TRUE,go=TRUE,nt=1)
## Gradient only BFGS optimizer to estimate smoothing parameters of models fitted by
## gam.fit1. The point of gradient only is to avoid evaluating log determinants during
## optimization. Set up for minimization!
## * stab (T/F) controls whether or not stabilizing transformations are to be applied.
## * go (T/F) controls whether to use gradiant only or regular BFGS.
## * L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
##   parameters actually multiplying the S[[i]]'s. sp's do not include the 
##   log scale parameter here.
##
## BFGS is based on Nocedal & Wright (2006) Numerical Optimization, Springer.
## In particular the step lengths are chosen to meet modified Wolfe conditions
## using their algorithms 3.5 (p60) and 3.6 (p61). The modification replaces 
## Wolfe 1 with a gradient only based alternative. 
## On p143 N&W recommend a post step
## adjustment to the initial Hessian. I can't understand why one would do anything
## other than adjust so that the initial Hessian would give the step taken, and
## indeed the latter adjustment seems to give faster convergence than their 
## proposal, and is therefore implemented.
##
{ zoom <- function(lo,hi) {
  ## local function implementing Algorithm 3.6 of Nocedal & Wright
  ## (2006, p61) Numerical Optimization. Relies on R scoping rules. 
  ## alpha.lo and alpha.hi are the bracketing step lengths.
  ## This routine bisection searches for a step length that meets the
  ## Wolfe conditions. lo and hi are both objects containing fields
  ## `score', `alpha', `dscore', where `dscore' is the derivative of 
  ## the score in the current step direction, `grad' and `mustart'. 
  ## `dscore' will be NULL if the gradiant has yet to be evaluated.
    for (i in 1:40) {
      trial <- list(alpha = (lo$alpha+hi$alpha)/2)
      lsp <- ilsp + step * trial$alpha
      b <- gam.fit1(G,Sl,sp=L%*%lsp+lsp0,scale,Sld=Sld,start=lo$start,mustart=lo$mustart,control = control,stab=stab,nt=nt)

      trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0)
      trial$grad <- t(L)%*%b$laml1;
      trial$mustart <- b$fitted.values
      trial$scale.est <- b$scale.est ## previously dev, but this differs from newton
      trial$start <- b$coefficients
      trial$score <- b$laml; 
      trial$dscore <- sum(step*trial$grad) ## directional derivative
      rm(b)
      notWolfe1 <- if (go) trial$dscore>-c1*initial$dscore||trial$dscore> c1*abs(lo$dscore) else
                  trial$score>initial$score+trial$alpha*c1*initial$dscore||trial$score>=lo$score 
      if (notWolfe1) {
        hi <- trial ## failed Wolfe 1
      } else { ## met Wolfe 1
        if (abs(trial$dscore) <= -c2*initial$dscore) return(trial) ## met Wolfe 2
        ## failed Wolfe 2 ...
        if (trial$dscore*(hi$alpha-lo$alpha)>=0) {
          hi <- lo }  
        lo <- trial 
      }  
    } ## end while(TRUE)
    return(NULL) ## failed
  } ## end zoom

  ## sanity check L
  if (is.null(L)) L <- diag(length(lsp)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(lsp0)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,nrow(L))
 
  ## initial fit...

  initial.lsp <- ilsp <- lsp
  b <- gam.fit1(G,Sl,sp=L%*%ilsp+lsp0,scale,Sld=Sld,start=start,mustart=mustart,control = control,stab=stab,nt=nt)

  initial <- list(alpha = 0,mustart=b$fitted.values,start=coef(b))
  score <- b$laml; grad <- t(L)%*%b$laml1;

  ## NOTE: following commented out does not match gam.fit1 dVkk returned for all...
  ## dVkk only refers to smoothing parameters, but sp may contain
  ## extra parameters at start and scale parameter at end. Have
  ## to reduce L accordingly... 
  #if (!is.null(G$family$n.theta)&&G$family$n.theta>0) {
  #  ind <- 1:G$family$n.theta
  #  nind <- ncol(L) - G$family$n.theta - if (G$family$n.theta + nrow(b$dVkk)<nrow(L)) 1 else 0 
  #  spind <- if (nind>0) G$family$n.theta+1:nind else rep(0,0)
  #  rspind <- G$family$n.theta + 1:nrow(b$dVkk)
  #} else {
    nind <- ncol(L) - if (nrow(b$dVkk)<nrow(L)) 1 else 0 
    spind <- if (nind>0) 1:nind else rep(0,0) ## index of smooth parameters
    rspind <- 1:nrow(b$dVkk)
  #}  
  L0 <- L[rspind,spind] ##if (nrow(L)!=nrow(b$dVkk)) L[spind,spind] else L
  
  initial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0)
  initial$score <- score;initial$grad <- grad;
  initial$scale.est <- b$scale.est
  start0 <- coef(b)
  mustart0 <- fitted(b)
  rm(b)

  B <- diag(length(initial$grad)) ## initial Hessian
  feps <- 1e-4
  for (i in 1:length(lsp)) { ## loop to FD for Hessian
     ilsp <- lsp;ilsp[i] <- ilsp[i] + feps
     b <- gam.fit1(G,Sl,sp=L%*%ilsp+lsp0,scale,Sld=Sld,start=start0,mustart=mustart0,control = control,stab=stab,nt=nt)
     grad1 <- t(L)%*%b$laml1;
     B[i,] <- (grad1-grad)/feps 
     rm(b)
  } ## end of FD Hessian loop
  ## force initial Hessian to +ve def and invert... 
  B <- (B+t(B))/2
  eb <- eigen(B,symmetric=TRUE)
  eb$values <- abs(eb$values)
  thresh <- max(eb$values) * 1e-4
  eb$values[eb$values<thresh] <- thresh
  B <- eb$vectors%*%(t(eb$vectors)/eb$values)
  ilsp <- lsp
  max.step <- 200

  c1 <- 1e-4;c2 <- .9 ## Wolfe condition constants

  score.hist <- rep(NA,max.step+1)
  score.hist[1] <- initial$score

  check.derivs <- FALSE;eps <- 1e-5

  uconv.ind <- rep(TRUE,ncol(B))
  rolled.back <- FALSE 

  for (i in 1:max.step) { ## the main BFGS loop
   
    ## get the trial step ...
    step <- initial$grad*0
    step[uconv.ind] <- -B[uconv.ind,uconv.ind]%*%initial$grad[uconv.ind]

    ## following tends to have lower directional grad than above (or full version commented out below)
    #step <- -drop(B%*%initial$grad)
    ## following line would mess up conditions under which Wolfe guarantees update,
    ## *if* based only on grad and not grad and hess...  
    #step[!uconv.ind] <- 0 ## don't move if apparently converged 
    
    if (sum(step*initial$grad)>=0) { ## step not descending!
      ## Following would really be in the positive definite space... 
      ##step[uconv.ind] <- -solve(chol2inv(chol(B))[uconv.ind,uconv.ind],initial$grad[uconv.ind])
      step <- -diag(B)*initial$grad ## simple scaled steepest descent 
      step[!uconv.ind] <- 0 ## don't move if apparently converged 
    }

    ms <- max(abs(step))
    trial <- list()
    if (ms>maxNstep) { 
      trial$alpha <- maxNstep/ms
      alpha.max <- trial$alpha*1.05
      ## step <- maxNstep * step/ms
      #alpha.max <- 1 ## was 50 in place of 1 here and below
    } else {
      trial$alpha <- 1 
      alpha.max <- min(2,maxNstep/ms) ## 1*maxNstep/ms
    }
    initial$dscore <- sum(step*initial$grad)
    prev <- initial
 
    while(TRUE) { ## step length control Alg 3.5 of N&W (2006, p60)
      lsp <- ilsp + trial$alpha*step
      b <- gam.fit1(G,Sl,sp=L%*%lsp+lsp0,scale,Sld=Sld,start=prev$start,mustart=prev$mustart,control = control,stab=stab,nt=nt)
     
      ### Derivative testing code. Not usually called and not part of BFGS...
      ok <- check.derivs
      while (ok) { ## derivative testing
        ## WARNING: if you zero log|S| or log|H| for testing in gam.fit1 then turn of
	##          stabilization or derivatives can't match (terms have cancelling repara corrections)
        #deriv <- 1
        ok <- FALSE ## set to TRUE to re-run (e.g. with different eps)
       	 
        bb <- b 

        fd.laml1 <- bb$laml1
        fdb <- bb$db*0
        for (j in 1:length(lsp)) { ## check dH and db.drho
          lsp1 <- lsp;lsp1[j] <- lsp[j] + eps
          ba <- gam.fit1(G,Sl,sp=L%*%lsp1+lsp0,scale,Sld=Sld,start=prev$start,mustart=prev$mustart,control = control,stab=stab,nt=nt)
	 
          fd.laml1[j] <- (ba$laml - bb$laml)/eps
          if (j<=ncol(fdb)) fdb[,j] <- (ba$coefficients - bb$coefficients)/eps
        }
	plot(fdb,bb$db);abline(0,1)
      } 
      ### end of derivative testing. BFGS code resumes...
      trial$score <- b$laml;
      trial$grad <- t(L)%*%b$laml1;
      trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0) ## curvature testing matrix
      trial$dscore <- sum(trial$grad*step)
      trial$mustart <- b$fitted.values
      trial$start <- b$coefficients
      trial$scale.est <- b$scale.est
      
      rm(b)
      Wolfe2 <- TRUE
      ## check the first Wolfe condition (sufficient decrease)...
      notWolfe1 <- if (go) trial$dscore>-c1*initial$dscore||trial$dscore> c1*abs(prev$dscore) else
                   (trial$score>initial$score+c1*trial$alpha*initial$dscore)#||(deriv==0&&trial$score>=prev$score)
      if (notWolfe1) {
         trial <- zoom(prev,trial) ## Wolfe 1 not met so backtracking
         break
      } 
      
      ## Note that written this way so that we can pass on to next test when appropriate...
     
      if (abs(trial$dscore) <= -c2*initial$dscore) break; ## `trial' is ok. (2nd Wolfe condition met).
      Wolfe2 <- FALSE

      if (!go && trial$dscore>=0) { ## increase at end of trial step
        trial <- zoom(trial,prev)
        Wolfe2 <- if (is.null(trial)) FALSE else TRUE
        break
      }
      
      prev <- trial
      if (trial$alpha == alpha.max) break ## { trial <- NULL;break;} ## step failed
      trial <- list(alpha = min(prev$alpha*1.3, alpha.max)) ## increase trial step to try to meet Wolfe 2
    } ## end of while(TRUE)

    ## Now `trial' contains a suitable step, or is NULL on complete failure to meet Wolfe,
    ## or contains a step that fails to meet Wolfe2, so that B can not be updated  
    if (is.null(trial)) { ## step failed
      lsp <- ilsp
      break ## failed to move, so nothing more can be done. 
    } else { ## update the Hessian etc...
     
      yg <- trial$grad-initial$grad
      step <- step*trial$alpha
      rho <- sum(yg*step)
      if (rho>0) { #Wolfe2) { ## only update if Wolfe2 is met, otherwise B can fail to be +ve def.
        if (i==1) { ## initial step --- adjust Hessian as p143 of N&W
          B <- B * trial$alpha ## this is my version 
          ## B <- B * sum(yg*step)/sum(yg*yg) ## this is N&W
        }
        rho <- 1/rho ## 1/sum(yg*step)
        B <- B - rho*step%*%(t(yg)%*%B)

        ## Note that Wolfe 2 guarantees that rho>0 and updated B is 
        ## +ve definite (left as an exercise for the reader)...
        B <- B - rho*(B%*%yg)%*%t(step) + rho*step%*%t(step)
      }

      score.hist[i+1] <- trial$score

      lsp <- ilsp <- ilsp + step 

      ## test for convergence
      converged <- TRUE
      score.scale <- 1 + abs(trial$score) ## abs(log(trial$dev/nrow(X))) + abs(trial$score) 
      uconv.ind <- abs(trial$grad) > score.scale*conv.tol 
      if (sum(uconv.ind)) converged <- FALSE
      ## following must be tighter than convergence...
      uconv.ind <- abs(trial$grad) > score.scale*conv.tol*.1 
      uconv.ind[spind] <- uconv.ind[spind] | abs(trial$dVkk) > score.scale * conv.tol*.1 
      if (abs(initial$score-trial$score) > score.scale*conv.tol) { 
        if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
        converged <- FALSE      
      }

      ## roll back any `infinite' smoothing parameters to the point at
      ## which score carries some information about them and continue 
      ## optimization. Guards against early long steps missing shallow minimum. 
      if (converged) { ## try roll back for `working inf' sps...
        if (sum(!uconv.ind)==0||rolled.back) break
        rolled.back <- TRUE
        counter <- 0
        uconv.ind0 <- uconv.ind 
        while (sum(!uconv.ind0)>0&&counter<5) {
          ## shrink towards initial values...
          lsp[!uconv.ind0] <- lsp[!uconv.ind0]*.8 + initial.lsp[!uconv.ind0]*.2
	  b <- gam.fit1(G,Sl,sp=L%*%lsp+lsp0,scale,Sld=Sld,start=trial$start,mustart=trial$mustart,control = control,stab=stab,nt=nt)
          trial$score <- b$laml
          trial$grad <- t(L)%*%b$laml1;
          trial$dscore <- sum(trial$grad*step)
          trial$scale.est <- b$scale.est
          trial$dVkk <- diag(t(L0) %*% b$dVkk %*% L0) ## curvature testing matrix 
          rm(b);counter <- counter + 1
          ## note that following rolls back until there is clear signal in derivs...
          uconv.ind0 <- abs(trial$grad) > score.scale*conv.tol*20        
          uconv.ind0[spind] <- uconv.ind0[spind] |  abs(trial$dVkk) > score.scale * conv.tol * 20
          uconv.ind0 <- uconv.ind0 | uconv.ind ## make sure we don't start rolling back unproblematic sps 
        }
        uconv.ind <- uconv.ind | TRUE
        ## following line is tempting, but will likely reduce usefullness of B as approximation 
        ## to inverse Hessian on return...
        ##B <- diag(diag(B),nrow=nrow(B))
        ilsp <- lsp
      }
    
      initial <- trial
      initial$alpha <- 0
    }  
  } ## end of iteration loop


  if (is.null(trial)) { 
    ct <- "step failed"
    lsp <- ilsp
    trial <- initial
  }
  else if (i==max.step) ct <- "iteration limit reached" 
  else ct <- "full convergence"
  ## final fit
  b <- gam.fit1(G,Sl,sp=L%*%lsp+lsp0,scale,Sld=Sld,start=trial$start,mustart=trial$mustart,control = control,stab=stab,nt=nt)
  score <- b$laml;grad <- t(L)%*%b$laml1;
 
  b$dVkk <- NULL
  ## get approximate Hessian (why???)...
  #ev <- eigen(B,symmetric=TRUE)
  #ind <- ev$values>max(ev$values)*.Machine$double.eps^.9
  #ev$values[ind] <- 1/ev$values[ind]
  #ev$values[!ind] <- 0
  #Bi <- ev$vectors %*% (ev$values*t(ev$vectors))

  list(score=score,lsp=lsp,lsp.full=L%*%lsp+lsp0,grad=grad,ihess=B,iter=i,conv =ct,
       score.hist=score.hist[!is.na(score.hist)],object=b)
} ## end of go.bfgs



inisp <- function(G,Sl,etastart=NULL,mustart=NULL,start=NULL,nt=1) {
## initial smoothing parameters to be found before any reparameterization
## of G$X... Operating post reparamaterization would not allow use of initial.sp
  if (length(G$S)==0) return(rep(0,0))
  family <- G$family
  y <- G$y
  weights <- G$w
  nobs <- if (is.matrix(y)) nrow(y) else length(y)
  nop <- ncol(G$X)
  discrete <-  !is.null(G$Xd)
  x <- if (discrete) list(Xd = G$Xd,lpid=G$lpid,kd = G$kd,ks=G$ks,ts=G$ts,dt=G$dt,qc=G$qc,v=G$v,p=nop) else G$X
  if (discrete) {
    attr(x,"lpi") <- attr(G$X,"lpi")
    if (inherits(G$family,"general.family")&&is.null(G$family$discrete.ok)) stop("family does not support discrete computation methods")
  }
  E <- attr(Sl,"E") ## E'E is a regularizing penalty which may be used in initialization
  ## E needs to be back in original parameterization here
  E <- Sl.initial.repara(Sl,E,inverse=TRUE,both.sides=FALSE,cov=FALSE)
  attr(E,"use.unscaled") <- TRUE 
  offset <- G$offset
  start0 <- start;start <- NULL	
  if (is.null(mustart)) {
        eval(family$initialize)
  } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
  }
  
  if (!is.null(start0)) start <- start0
  
  if (inherits(family,"general.family")) { ## Cox, gamlss etc...   
    lbb <- -family$ll(y,x,start,G$weights,family,offset=G$offset,deriv=1)$lbb
    theta <- rep(0,0)
  } else { ## extended or exponential families...
    if (is.null(mustart)) {
      if (!is.null(start)) etastart <- offset + if (discrete) mgcv:::Xbd(G$Xd,start,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) else drop(x%*%start)
      if (!is.null(etastart)) mustart <- family$linkinv(etastart)
    } #else mustart <- mukeep
    if (inherits(family,"extended.family")) {
      theta <- family$getTheta()
      Ddo <- family$Dd(y,mustart,theta,G$w)
      mu.eta2 <-family$mu.eta(family$linkfun(mustart))^2 
      w <- .5 * as.numeric(Ddo$Dmu2 * mu.eta2)
    } else w <- as.numeric(G$w*family$mu.eta(family$linkfun(mustart))^2/family$variance(mustart))
    w[!is.finite(w)] <- 0
    lbb <- if (discrete) mgcv:::XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,nt,G$drop) else crossprod(x,w*x)
  }  
  lsp <- log(abs(initial.sp(lbb,G$S,G$off,XX=TRUE)))
  ## Deal with L, scale and theta ...
  if (!is.null(G$L)) lsp <- if (ncol(G$L)>0) as.numeric(coef(lm(lsp ~ G$L-1+offset(G$lsp0)))) else rep(0,0)
  if (inherits(family,"extended.family")) lsp <- c(theta,lsp)
  if (G$scale<=0) lsp <- c(lsp,log(var(y)*.05))
  ## Following is a *really* stupid idea - need not be +ve definite at this point.
  ## now create penalized lbb suitable for detecting rank deficiency...
  #lbb <- lbb + crossprod(E) ## penalized in original space
  #lbb <- Sl.initial.repara(Sl,lbb,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=nt) ## convert to initially reparameterized
  ## now run a pivoted Cholesky to actually perform the check
  
  lsp
} ## inisp

go.estimate <- function(G,stab=TRUE,go=FALSE,etastart=NULL,mustart=NULL,start=NULL,scale = -1) {
## Wrapper to estimate gam via go.bfgs/gam.fit1 combination
## * G is a pre-fit gam object
## * stab (T/F) controls whether to use stability transforms.
## * go (T/F) control whether to use gradient only or regular BFGS
## stab and go both FALSE is not recommended
  n.theta <- 0
  if (inherits(G$family,"extended.family")) if (inherits(G$family,"general.family")) {
    if (!is.null(G$family$preinitialize)) {
      Gmod <- G$family$preinitialize(G) ## modifies some elements of G
      #for (gnam in names(Gmod)) G[[gnam]] <- Gmod[[gnam]] ## copy these into G
      G[names(Gmod)] <- Gmod
    }
    scale <- 1
  } else { ## preinitialize extended family
      pini <- if (is.null(G$family$preinitialize)) NULL else G$family$preinitialize(y,G$family)
      if (!is.null(pini$Theta)) G$family$putTheta(pini$Theta)
      if (!is.null(pini$y)) y <- pini$y
      if (scale==0) scale <- -1
      if (is.null(G$family$scale)) scale <- 1 else scale <- if (G$family$scale<0) scale else G$family$scale
      n.theta <- G$family$n.theta
  }
  if (G$family$family %in% c("binomial","poisson") && scale <= 0) scale <- 1
  G$scale <- scale
  Sl <- Sl.setup(G,cholesky=TRUE) ## set up efficient smooth list
  lsp <- inisp(G,Sl,etastart=etastart,mustart=mustart,start=start) ## get the initial smoothing parameters
  discrete <- !is.null(G$Xd)
  if (!discrete) G$X <- Sl.initial.repara(Sl,G$X,inverse=FALSE,both.sides=FALSE,cov=FALSE,nt=1)
  G$family <- fix.family.ls(fix.family.link(fix.family.var(G$family)))
 
  L <- if (is.null(G$L)) diag(length(G$lsp0)) else G$L
  lsp0 <- G$lsp0
  if (n.theta) {
    lsp0 <- c(rep(0,n.theta),lsp0)
    L <- cbind(rbind(diag(n.theta),matrix(0,nrow(L),n.theta)),rbind(matrix(0,n.theta,ncol(L)),L))
  }
  if (scale<0) {
    lsp0 <- c(lsp0,0)
    L <- rbind(cbind(L,0),0)
    L[nrow(L),ncol(L)] <- 1
  }

  ## at this stage it is necessary to deal with structural lack of identifiability...
  ## first run a fit and establish structural identifiability without stability reparameterization
  ## idea is that at the initial sps any rank deficiency had better be structural or model is
  ## deeply wrong, but we want to drop in the space of the initial reparameterization...
  
  idrop <- gam.fit1(G,Sl,sp=drop(L%*%lsp+lsp0),scale=G$scale,stab=FALSE,ident=TRUE)
  np <- ncol(G$X)
  if (length(idrop)) {
    Sld <- mgcv:::Sl.drop(Sl,idrop,np) ## note idrop added as 'drop' attribute of Sld
    if (!discrete) {
      lpi <- attr(G$X,"lpi")
      if (!is.null(lpi)) { ## needs to be modified to allow for explicit col dropping
        idrop <- sort(idrop,decreasing=TRUE)
        for (i in idrop) for (j in 1:length(lpi)) { ## adjust the indices for dropped terms
          lpi[[j]] <- lpi[[j]][lpi[[j]]!=i]
          lpi[[j]][lpi[[j]]>i] <- lpi[[j]][lpi[[j]]>i] - 1
        }
      }	
      G$X <- G$X[,-idrop]
      attr(G$X,"lpi") <- lpi
    }
    if (!is.null(G$g.index)) {
      if (any(G$g.index[idrop])) stop("non-linear parameters appear unidentifiable - no way to proceed")
      G$g.index <- G$g.index[-idrop]
    }
  } else Sld <- NULL
  ## Sld never used for initial repara - Sl.drop does not modify the repara records
  
## NOTE: I think any supplied start must need to be reparameterized 

  chder <- FALSE
  while (chder) { ## debugging derivative testing code (assumes L = I)
    chder <- FALSE
    eps <- 1e-4
    er <- gam.fit1(G,Sl,sp=lsp,scale=G$scale,Sld=Sld,stab=stab)
    fd.laml <- er$laml1*0
    fd.db <- er$db*0
    fd.ddetH <- er$ddetH*0
    for (i in 1:length(lsp)) {
      lsp1 <- lsp;lsp1[i] <- lsp[i] + eps
      um <- gam.fit1(G,Sl,sp=lsp1,scale=G$scale,Sld=Sld,stab=stab)
      fd.laml[i] <- (um$laml - er$laml)/eps
      if (i<=ncol(fd.db)) {
        fd.db[,i] <- (um$beta - er$beta)/eps
        fd.ddetH[i] <- (um$detH - er$detH)/eps
        cat(".")
        plot(er$db[,i],fd.db[,i],pch=19,cex=.3,main=i);abline(0,1)
      }
    }
    print(rbind(er$laml1,fd.laml))
    plot(er$db,fd.db,pch=19,cex=.3);abline(0,1)
  } ## chder

 
  um <- go.bfgs(lsp=lsp,G,Sl,scale=G$scale,lsp0=lsp0,Sld=Sld,L=L,start=NULL,etastart=NULL,
                     mustart=NULL,control = gam.control(),stab=stab,go=go,nt=1)

  attr(Sl,"drop") <- idrop
  ud <- gam.fit1.postproc(um$object,Sl,S=G$S,off=G$off,V.rho=um$ihess,L=L,lsp=um$lsp)
  #um[names(ud)] <- ud ## lovely!
  ## Create a gam object: (scale is Fletcher estimate)
  names(G$cmX) <- names(ud$edf2) <- names(ud$edf1) <- names(ud$edf) <- names(ud$coef) <- G$term.names
  object <- list(coefficients=ud$coef,family=G$family,y=G$y,prior.weights=G$w,fitted.values=um$object$fitted.values,
            linear.predictors=um$object$fitted.values,rank=ud$rank,aic=um$object$aic+sum(ud$edf),
            gcv.ubre=um$score,sp=exp(lsp),outer.info=list(conv=um$conv,iter=um$iter,grad=um$grad,score.hist=um$score.hist),
            scale.estimated=!um$object$scale.known,scale=ud$scale.est,method="REML",Vp=ud$Vb,Vc=ud$Vc,Ve=ud$Ve,edf=ud$edf,
	    edf1=ud$edf1,edf2=ud$edf2,R=ud$R,nsdf=G$nsdf,sig2=um$object$scale.est,smooth=G$smooth, full.sp=exp(um$lsp.full),
	    formula=G$formula,var.summary=G$var,summary,cmX=G$cmX,model=G$mf,terms=G$terms,pred.formula=G$pred.formula,
            pterms=G$pterms,assign=G$assign,xlevels=G$xlevels,offset=G$offset,df.residual=G$n-sum(ud$edf),
            min.edf=G$min.edf,optimizer="go.bfgs",call=G$cl,deviance=um$object$deviance)

  if (!inherits(G$family,"extended.family")) object$null.deviance <-
         glm(object$y~offset(G$offset),family=object$family,weights=object$prior.weights)$deviance
	 
  ## extended family may need to manipulate fit object. Code
  ## will need to include the following line if G$X used...
  ## G$X <- Sl.initial.repara(G$Sl,G$X,inverse=TRUE,cov=FALSE,both.sides=FALSE)
  if (!is.null(G$family$postproc)) {
    if (inherits(G$family,"general.family")) eval(G$family$postproc) else {
      posr <- G$family$postproc(family=object$family,y=G$y,prior.weights=object$prior.weights,
              fitted=object$fitted.values,linear.predictors=object$linear.predictors,offset=G$offset,
	      intercept=G$intercept)	      
      if (!is.null(posr$family)) object$family$family <- posr$family
      if (!is.null(posr$deviance)) object$deviance <- posr$deviance
      if (!is.null(posr$null.deviance)) object$null.deviance <- posr$null.deviance	      
    }
  }
  
  if (is.list(object$formula)) attr(object$formula,"lpi") <- attr(G$X,"lpi")
  attr(object$pred.formula,"full") <- reformulate(all.vars(object$terms))

  ## in the discrete case, save the discretization information for use in efficient prediction...
  if (discrete) object$dinfo <- list(gp=G$gp, v = G$v, ts = G$ts, dt = G$dt, qc = G$qc, drop = G$drop,
                                     lpid=G$lpid,lpip=attr(G$Xd,"lpip"))

  ## avoid large environments being stored with object.... 
  environment(object$formula) <- environment(object$pred.formula) <-
  environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
  if (!is.null(object$model))  environment(attr(object$model,"terms"))  <- .GlobalEnv
  if (!is.null(attr(object$pred.formula,"full"))) environment(attr(object$pred.formula,"full")) <- .GlobalEnv

  class(object) <- c("gam","glm","lm")
  if (is.null(object$deviance)) object$deviance <- sum(residuals(object,"deviance")^2)

  object
} ## go.estimate


smooth.construct.up.smooth.spec <- function(object,data,knots) {
  object$mono = 1;
  smooth.construct.ps.smooth.spec(object,data,knots)
}