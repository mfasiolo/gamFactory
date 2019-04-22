#' Creates stack effects for a SHASH family
#'
#' @description Creates a family to be used in mgcv::gam to estimate a generalized additive model
#' with response variable having SHASH distribution and stack nonnegative effects. 
#' 
#' @param X Model Matrix of the nonnegative linear combination
#' @param link At the moment only the identity link is possible
#'
#' @return A family to give as family argument to mgcv::gam
#' @export
#'
#' @examples
stackFamilySHASH <- function(X, link="identity") {
  
  K <- ncol(X)
  link <- lapply(1:K, function(x) "identity")
  stats <- list()
  for (ii in 1:K) {
    stats[[ii]] <- make.link(link[[ii]])
    fam <- structure(list(link=link[[ii]],canonical="none",linkfun=stats[[ii]]$linkfun,
                          mu.eta=stats[[ii]]$mu.eta), class="family")
    fam <- fix.family.link(fam)
    stats[[ii]]$d2link <- fam$d2link
    stats[[ii]]$d3link <- fam$d3link
    stats[[ii]]$d4link <- fam$d4link
  }
  
  residuals <- function(object, type) {
    return(y)
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
    ntheta <- 3 ## number of additional parameters of the SHASH distribution
    lpi <- attr(G$X,"lpi")
    #offs <- attr(G$X,"offset")
    XX <- crossprod(G$X)
    G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta)) ## add dummy columns to G$X
    #G$cmX <- c(G$cmX,rep(0,ntheta)) ## and corresponding column means
    G$term.names <- c(G$term.names,paste("R",1:ntheta,sep="."))
    attr(G$X,"lpi") <- lpi
    
    
    #offs -> attr(G$X,"offset")
    attr(G$X,"XX") <- XX
    
    # ## pad out sqrt of balanced penalty matrix to account for extra params
    
    ### !!!!!!!!!!!!!!! This is not returned by the function, maybe because before it was an expression???
    # if (!is.null(G$Sl)) attr(G$Sl,"E") <- cbind(attr(G$Sl,"E"),matrix(0,nbeta,ntheta)) 
    
    G$family$data <- list(ydim = ydim,nbeta=nbeta) # !!!!!!!!!!!!!!!!!!!!!!!!! Do you need this???
    G$family$ibeta = rep(0,ncol(G$X))
    # ## now get initial parameters and store in family...
    # for (k in 1:ydim) {
    #   sin <- G$off %in% lpi[[k]]
    #   #Sk <- G$S[sin]
    #   um <- magic(G$y[,k],G$X[,lpi[[k]]],rep(-1,sum(sin)),G$S[sin],
    #               match(G$off[sin],lpi[[k]])) # , ## turn G$off global indices into indices for this predictor
    #   #nt=control$nthreads)
    #   G$family$ibeta[lpi[[k]]] <- um$b
    #   G$family$ibeta[nbeta+1] <- -.5*log(um$scale) ## initial log root precision
    #   nbeta <- nbeta + ydim - k + 1
    # }
    
    
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
    K <- ncol(X)
    nu <- list()
    for (kk in 1:K) {
      nu[[kk]] <- Z %*% coef[lpi[[kk]]]
    }
    nu <- do.call("cbind", nu)
    a <- exp(nu)
    eta <- rowSums(a * X)
    mu <- eta
    theta.index <- (length(coef) - 2):length(coef)
    theta <- coef[theta.index]
    beta <- coef[- theta.index]
    pars <- cbind(mu, matrix(rep(theta, nn), nrow = nn, byrow = T))
    objSH <- createSH(y = y)$derObj(param = pars, deriv = deriv + 1)
    
    # log-likelihood
    l <- objSH$d0(SUM = T)
    
    if (deriv > 0) {
      # 1st derivatives
      le <- objSH$d1(SUM = F)[[1]]   # nx1 1st deriv wrt eta
      lt <- objSH$d1(SUM = T)[][- 1] # 3x1 1st deriv wrt theta
      
      # 2nd derivatives
      lee <- objSH$d2(SUM = F)[[1]]  # nx1 2nd deriv wrt eta
      let <- do.call(cbind, objSH$d2(SUM = F)[2:4]) # nx3 2nd mixed deriv wrt eta theta
      ltt <- matrix(NA, 3, 3)          # 3x3 2nd deriv wrt theta
      ltt[lower.tri(ltt, diag=TRUE)] <- objSH$d2(SUM = T)[- (1:4)]
      ltt <- pmax(ltt, t(ltt), na.rm=TRUE)
    }
    
    if (deriv > 1) {
      
      #3rd derivatives
      leee <- objSH$d3(SUM = F)[[1]] # nx1 3rd deriv wrt eta
      leet <- do.call(cbind, objSH$d3(SUM = F)[2:4]) # nx3 3rd deriv wrt eta eta theta
      lett <- do.call(cbind, objSH$d3(SUM = F)[5:10]) # nx(3*(3+1)/2) 3rd deriv wrt eta theta theta
      lttt <- objSH$d3(SUM = T)[11:20] # array of 3rd deriv wrt theta
      
      i3 <- trind.generator(3)$i3
      ltttVec <- objSH$d3(SUM = T)[11:20] # 3x3x3 array of 3rd deriv wrt theta
      lttt <- array(NA, dim = c(3, 3, 3))
      for (uu in 1:3) for (vv in 1:3) for (ww in 1:3) {
        lttt[uu, vv, ww] <- ltttVec[i3[uu, vv, ww]]
        ww <- ww + 1
      }
      
    }
    
    if (deriv) {
      ## TBD: transform derivates w.r.t. mu to derivatives w.r.t. eta...
      
      ret <- convertDerivStackPositive(beta, X, Z, 
                                       le, lt, 
                                       lee, let, ltt, 
                                       leee = leee, leet = leet, 
                                       lett = lett, lttt = lttt, 
                                       d1b = d1b, deriv = deriv - 1)
    } else ret <- list()
    ret$l <- l
    ret
  } # end ll stackFamilySHASH
  
  

  
  rd <- function(mu,wt,scale) {
    
  } ## rd
  
  dev.resids <- function(a, b, c, d) y # MAYBE IT'S NEEDED IN gam.fit5
  
  structure(list(family="stackSHASH",ll=ll,nlp=K,
                 link="identity",
                 preinitialize=preinitialize,
                 initialize=initialize,
                 # postproc=postproc,residuals=residuals,
                 linfo = stats,
                 rd=rd,
                 dev.resids = dev.resids,
                 linkinv = stats$linkinv, # MAYBE IT'S NEEDED IN gam.fit5
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 1 ## signal only first derivatives available...
  ), class = c("general.family","extended.family","family"))
} # stackFamilySHASH

