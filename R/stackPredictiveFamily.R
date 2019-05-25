#' Creates a family for stacking predictive distributions
#'
#' @description Creates a family to be used in mgcv::gam 
#' for stacking predictive distributions.
#' 
#' @param logP Matrix of the log-predictive distributions to be stacked
#'
#' @return A family to give as family argument to mgcv::gam
#' @export
#' 
#' @name stackPredictiveFamily
#' @rdname stackPredictiveFamily
#'
stackPredictiveFamily <- function(logP) {
  
  link <- "identity"
  K <- ncol(logP)
  link <- lapply(1:(K - 1), function(x) "identity")
  stats <- list()
  for (ii in 1:(K - 1)) {
    stats[[ii]] <- make.link(link[[ii]])
    fam <- structure(list(link=link[[ii]],canonical="none",linkfun=stats[[ii]]$linkfun,
                          mu.eta=stats[[ii]]$mu.eta), class="family")
    fam <- fix.family.link(fam)
    stats[[ii]]$d2link <- fam$d2link
    stats[[ii]]$d3link <- fam$d3link
    stats[[ii]]$d4link <- fam$d4link
  }
  
  residuals <- function(object, type=c("deviance","pearson","response")) {
    
    return(object$y[, 1])
    # type <- match.arg(type)
    # if (type %in% c("deviance", "pearson")) return(NULL)
    # else object$y[, 1] - rowSums(exp(object$fitted) * X)
  }  
  # residuals <- function(object,type=c("deviance","pearson","response")) {
  #   type <- match.arg(type)
  #   rsd <- object$y-object$fitted[,1]
  #   if (type=="response") return(rsd) else
  #     return((rsd*object$fitted[,2])) ## (y-mu)/sigma 
  # }
  
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## in principle the following seems reasonable, but because no
    ## price is paid for the high null variance, it leads to silly
    ## % deviance explained...
    #er <- fitNull(G$y,G$family,G$w,G$offset,nlp=length(attr(G$X,"lpi")),tol=1e-7)
    #object$null.deviance <- sum(((object$y-er$mu[,1])*er$mu[,2])^2*G$w)
    
    
    ### I HAVE NOT CHANGED ANYTHING YET HERE!!!
    object$null.deviance <- sum(((object$y-mean(object$y))*object$fitted[,2])^2)
  })
  
  
  preinitialize <- function(G) {
    ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
    ## elements, returning a named list of the modified ones.
    ## extends model matrix with dummy columns and 
    ## finds initial coefficients
    ydim <- ncol(G$y) ## dimension of response
    nbeta <- ncol(G$X)
    lpi <- attr(G$X,"lpi")
    #offs <- attr(G$X,"offset")
    XX <- crossprod(G$X)
    
    #offs -> attr(G$X,"offset")
    attr(G$X,"XX") <- XX
    
    # ## pad out sqrt of balanced penalty matrix to account for extra params
    
    G$family$data <- list(ydim = ydim, nbeta = nbeta)
    ## initialize all alphas equal by putting all betas to 0
    G$family$ibeta = rep(0,ncol(G$X)) 
    
    list(X=G$X,term.names=G$term.names,family=G$family)
  } ## preinitialize
  
  
  initialize <- expression({
    
    # n <- rep(1, nobs)
    if (is.null(start)) start <- family$ibeta
    
    ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WHAT happens here???
    # if (exists("rp", inherits = FALSE) && length(rp$rp) > 0) 
    #   attr(x, "XX") <- Sl.repara(rp$rp, t(Sl.repara(rp$rp, 
    #                                                 attr(x, "XX"))))
    
  }) ## initialize
  
  
  ll <- function(y,x,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    
    y <- y[, 1] # Using multiple linear predictors, y is repeated, here we do not want this.
    nn <- length(y)
    lpi <- attr(x,"lpi") ## extract linear predictor index, in gamlss.gH it's jj
    Z <- x[ , lpi[[1]], drop = FALSE]
    p <- ncol(Z); K <- ncol(logP)
    nu <- list()
    for (kk in 1:(K - 1)) nu[[kk]] <- Z %*% coef[lpi[[kk]]]
    nu <- do.call("cbind", nu)
    
    d1H <- lb <- lbb <- NULL ## default

    a <- cbind(1, exp(nu)) / (1 + rowSums(exp(nu)))
    am <- a[, - 1, drop = FALSE]
    
    ## adjust log-densities for numerical stability
    c <- apply(logP, 1, max)
    X <- logP - c
    Xm <- X[, - 1, drop = FALSE]
    w <- exp(nu + Xm - log(exp(X[, 1]) + rowSums(exp(nu + Xm))))
    
    ll <- log(rowSums(a * exp(logP)))

    if (deriv > 0) { ## grad and Hess
      ## the gradient...
      ln <- w - am
      lb <- t(Z) %*% ln
      lb <- as.vector(lb)
      
      ## the Hessian...
      lnn <- list()
      coun <- 0
      for (jj in 1:(K - 1)) for (kk in jj:(K - 1)) {
        coun <- coun + 1
        lnn[[coun]] <-
          ln[, jj] * (as.numeric(jj == kk) - w[, kk]) - am[, jj] * ln[, kk]
      }
      lnn <- do.call(cbind, lnn)
      i2 <- trind.generator(K - 1)$i2
      
      lbb <- lapply(1:(K - 1), function(x) vector("list", K - 1))
      for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
        lbb[[rr]][[ss]] <- t(Z * lnn[, i2[rr, ss]]) %*% Z
        if (ss > rr) lbb[[ss]][[rr]] <- lbb[[rr]][[ss]]
      }
      lbb <- lapply(lbb, function (x) do.call(rbind, x))
      lbb <- matrix(do.call(cbind, lbb), nrow = p * (K - 1))
    } ## grad and Hess
    
    if (deriv == 3) { ## store full d1H
      
      lnnn <- list()
      coun <- 0
      for (jj in 1:(K - 1)) for (kk in jj:(K - 1)) for (ll in kk:(K - 1)) {
        coun <- coun + 1
        lnnn[[coun]] <-
          lnn[, i2[jj, ll]] * (as.numeric(jj == kk) - w[, kk]) -
          w[, kk] * (as.numeric(kk == ll) - w[, ll]) * ln[, jj] -
          am[, jj] * (as.numeric(jj == ll) - am[, ll]) * ln[, kk] -
          am[, jj] * lnn[, i2[kk, ll]]
      }
      lnnn <- do.call("cbind", lnnn)
      
      m <- ncol(d1b)
      
      d1H <- list()
      for (l in 1:m) {
        d1brl <- matrix(d1b[, l], nrow = p, ncol = K - 1)
        lbbr <- lapply(1:(K - 1), function(x) vector("list", K - 1))
        i3 <- trind.generator(K - 1)$i3
        for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
          V <- rowSums((Z %*% d1brl) * lnnn[, i3[rr, ss, ]])
          lbbr[[rr]][[ss]] <- t(Z) %*% (Z * V)
          if (ss > rr) lbbr[[ss]][[rr]] <- lbbr[[rr]][[ss]]
        }
        lbbr <- lapply(lbbr, function (x) do.call(rbind, x))
        d1H[[l]] <- matrix(do.call(cbind, lbbr), 
                           nrow = p * (K - 1), 
                           ncol = p * (K - 1))
      } ## for (l in 1:m)
    } ## store full d1H
    
    list(l = sum(ll), lb = lb, lbb = lbb, d1H = d1H)
  } # end ll stackPredictiveFamily
  
  
  rd <- function(mu,wt,scale) {
    
  } ## rd
  
  dev.resids <- function(a, b, c, d) y # MAYBE IT'S NEEDED IN gam.fit5
  
  structure(list(family="stackPredictiveFamily",ll=ll,nlp=K - 1,
                 link="identity",
                 preinitialize=preinitialize,
                 initialize=initialize,
                 # postproc=postproc,
                 residuals=residuals,
                 linfo = stats,
                 rd=rd,
                 dev.resids = dev.resids,
                 linkinv = stats$linkinv, # MAYBE IT'S NEEDED IN gam.fit5
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 1 ## signal only first derivatives available...
  ), class = c("general.family","extended.family","family"))
} # stackFamily

