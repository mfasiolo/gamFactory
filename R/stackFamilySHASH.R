#' Title
#'
#' @param X 
#' @param link 
#'
#' @return
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
    object$null.deviance <- sum(((object$y-mean(object$y))*object$fitted[,2])^2)
  })
  
  ll <- function(y,x,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    
    y <- y[, 1] # Using multiple linear predictors, y is repeated, here we do not want this.
    jj <- attr(x,"lpi") ## extract linear predictor index
    K <- ncol(X)
    nu <- list()
    for (kk in 1:K) {
      Z <- x[ , jj[[kk]], drop = FALSE]
      nu[[kk]] <- Z %*% coef[jj[[kk]]]
    }
    nu <- do.call("cbind", nu)
    a <- exp(nu)
    eta <- rowSums(a * X)
    mu <- eta
    tau <- 1
    eps <- 1
    phi <- 1
    sig <- exp(tau)
    del <- exp(phi)
    theta <- c(tau, eps, phi)
    pars <- cbind(mu, theta[1], theta[2], theta[3])
    objSH <- createSH(y = y)$derObj(param = pars, deriv = 3)
    l <- objSH$d0(SUM = T)
    le <- objSH$d1(SUM = F)[[1]]
    lee <- objSH$d2(SUM = F)[[1]]
    leee <- objSH$d3(SUM = F)[[1]]
    ret <- convertDerivStackPositive(coef, theta, X, Z, le, lee, leee, d1b, 1)
    ret$l <- l
    ret
  }
  
  initialize <- expression({
    if (is.null(start)) {
      start <- rep(0, ncol(x))
    }
  }) ## initialize
  
  rd <- function(mu,wt,scale) {
    
  } ## rd
  
  dev.resids <- function(a, b, c, d) y # MAYBE IT'S NEEDED IN gam.fit5
  
  structure(list(family="stackSHASH",ll=ll,nlp=K,
                 link="identity",
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

