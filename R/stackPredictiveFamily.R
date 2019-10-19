#' Creates a family for stacking predictive distributions
#'
#' @description Creates a family to be used in mgcv::gam 
#' for stacking predictive distributions.
#' 
#' @param logP n x K matrix of the log-predictive distributions of the single experts. 
#' n is the total number of available observations, while K is the number of experts
#'
#' @return A family to give as family argument to mgcv::gam
#' @export
#' 
#' @name stackPredictiveFamily
#' @rdname stackPredictiveFamily
#' @examples 
#' library(gamFactory)
#' library(mgcViz)
#' 
#' # Generating training and stacking sets of increasing size
#' sizes <- c(rep(100, 60), rep(125, 20), rep(150, 5), 
#'            rep(175, 5), rep(200, 5), rep(300, 5),
#'            rep(400, 5), rep(600, 5), rep(800, 5))
#' 
#' ns <- length(sizes)
#' 
#' datListTrain <- lapply(sizes, 
#'                        function( .n ){
#'                          gamSim(1, n = .n, dist = "binary", 
#'                                 scale = 0.33, verbose = FALSE)
#'                        })
#' datTrain <- do.call("rbind", datListTrain)
#' datTrain <- as.data.frame(datTrain)
#' datTrain$sizes <- rep(sizes, sizes)
#' 
#' datListStack <- lapply(sizes, 
#'                        function( .n ){
#'                          gamSim(1, n = .n, dist = "binary", scale = 0.33, verbose = FALSE)
#'                        })
#' datStack <- do.call("rbind", datListStack)
#' datStack <- as.data.frame(datStack)
#' datStack$sizes <- rep(sizes, sizes)
#' 
#' # Estimate simple model and produce predictions on stacking set
#' m1 <- lapply(1:ns, 
#'              function(.kk){
#'                .fit <- gam(y ~ s(x0, k = 3) + s(x1, k = 3) + 
#'                                s(x2, k = 3) + s(x3, k = 3), 
#'                            family=binomial, data = datListTrain[[.kk]], 
#'                            method="REML")
#'                .y <- datListStack[[.kk]]$y
#'                .p <- predict(.fit, newdata = datListStack[[.kk]], type = "response")
#'                .stack <- log(.p) * .y + log1p( - .p ) * ( !.y )
#'                return( list("fit" = .fit, "stack" = .stack) )
#'              })
#' 
#' ### Estimate simple model and produce predictions on stacking set
#' m2 <- lapply(1:ns, 
#'              function(.kk){
#'                .fit <- gam(y ~ s(x0, k = 10) + s(x1, k = 10) +
#'                                s(x2, k = 10) + s(x3, k = 10), 
#'                            family=binomial, data = datListTrain[[.kk]], 
#'                            method="REML")
#'                .y <- datListStack[[.kk]]$y
#'                .p <- predict(.fit, newdata = datListStack[[.kk]], type = "response")
#'                .stack <- log(.p) * .y + log1p( - .p ) * ( !.y )
#'                return( list("fit" = .fit, "stack" = .stack) )
#'              })
#' 
#' # Perform the probabilistic stacking
#' logPStack <- cbind(do.call("c", lapply(m1, "[[", "stack")), 
#'                    do.call("c", lapply(m2, "[[", "stack")))
#' 
#' fitStack <- gamV(list(y ~ s(log(sizes), k = 7)), 
#'                  family = stackPredictiveFamily(logP = logPStack), 
#'                  data = datStack)
#' 
#' # The weight of the second model should increase with the size of the data set
#' plot(fitStack) + l_fitLine() + l_ciLine() + l_rug()
#' summary(fitStack)
#'
stackPredictiveFamily <- function(logP, rho = NULL) {
  
  link <- "identity"
  K <- ncol(logP)
  ## adjust log-densities for numerical stability
  X <- logP - matrixStats::rowMaxs(logP)
  expX <- exp(X)
  P <- exp(logP) ## predictive densities (not logs)
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
    
    return(as.matrix(object$y)[, 1])
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
    
    discrete <- is.list(x)
    lpi <- attr(x,"lpi") ## extract linear predictor index, in gamlss.gH it's jj
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    nu <- sapply(1:(K - 1), 
                 function(ii) {
                   if (discrete) {
                     Xbd(x$Xd,coef,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,
                         drop=x$drop,lt=x$lpid[[ii]])
                   } else {
                     x[, lpi[[ii]], drop = FALSE] %*% coef[lpi[[ii]]]
                   }
                 })
    
    
    d1H <- lb <- lbb <- NULL ## default
    
    nu1 <- cbind(0, nu)
    nuCen <- nu1 - rowMaxs(nu1)
    a <- exp(nuCen) / rowSums(exp(nuCen))
    
    w <- exp(nu + X[, - 1, drop = FALSE] - 
               log(expX[, 1] + 
                     rowSums(exp(nu + X[, - 1, drop = FALSE]))))
    
    logLik <- sum(log(rowSums(a * P)))
    if (!is.null(rho)) {
      logLik <- logLik + sum(t(t(log(a)) * (rho - 1)))
    }
    
    if (deriv > 0) { ## grad and Hess
      ## the gradient...
      ln <- ln_r <- w - a[, - 1, drop = FALSE]
      if (!is.null(rho)) {
        ln_r <- ln_r + t(t(- a[, - 1, drop = FALSE] * (sum(rho) - K)) + 
                           rho[- 1] - 1)
      }
      
      ## the Hessian...
      lnn <- lnn_rr <- list()
      coun <- 0
      for (jj in 1:(K - 1)) for (kk in jj:(K - 1)) {
        coun <- coun + 1
        lnn[[coun]] <- lnn_rr[[coun]] <- 
          ln[, jj] * (as.numeric(jj == kk) - w[, kk]) - 
          a[, jj + 1] * ln[, kk]
        if (!is.null(rho)) {
          lnn_rr[[coun]] <- 
            lnn_rr[[coun]] - a[, jj + 1] *
            (as.numeric(jj == kk) - a[, kk + 1]) *
            (sum(rho) - K)
        }
      }
      lnn <- do.call(cbind, lnn)
      lnn_rr <- do.call(cbind, lnn_rr)
    } ## grad and Hess
    
    lnnn <- l4 <- 0 ## defaults
    tri <- family$tri ## indices to facilitate access to earlier results
    
    if (deriv > 1) { ## store full d1H
      
      lnnn <- list()
      coun <- 0
      for (jj in 1:(K - 1)) for (kk in jj:(K - 1)) for (tt in kk:(K - 1)) {
        coun <- coun + 1
        lnnn[[coun]] <- 
          lnn[, tri$i2[jj, tt]] * (as.numeric(jj == kk) - w[, kk]) -
          w[, kk] * (as.numeric(kk == tt) - w[, tt]) * ln[, jj] -
          a[, jj + 1] * (as.numeric(jj == tt) - a[, tt + 1]) * ln[, kk] -
          a[, jj + 1] * lnn[, tri$i2[kk, tt]]
        if (!is.null(rho)) {
          lnnn[[coun]] <- lnnn[[coun]] -
            (
              (as.numeric(jj == kk) - a[, kk + 1]) * 
                (as.numeric(jj == tt) * a[, jj + 1] - a[, jj + 1] * a[, tt + 1]) - 
                (as.numeric(kk == tt) * a[, kk + 1] - a[, kk + 1] * a[, tt + 1]) * 
                a[, jj + 1]
            ) * (sum(rho) - K)
        }
      }
      lnnn <- do.call("cbind", lnnn)
    } 
    
    if (deriv) {
      ## get the gradient and Hessian...
      ret <- gamlss.gH(x,lpi,ln_r,lnn_rr,tri$i2,l3=lnnn,i3=tri$i3,l4=l4,i4=tri$i4,
                       d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l <- logLik; ret
    
  } # end ll stackPredictiveFamily
  
  
  rd <- function(mu,wt,scale) {
    
  } ## rd
  
  dev.resids <- function(a, b, c, d) y # MAYBE IT'S NEEDED IN gam.fit5
  
  structure(list(family="stackPredictiveFamily",ll=ll,nlp=K - 1,
                 link="identity",
                 getRho = getRho,
                 getLogP = getLogP,
                 putRho = putRho,
                 putLogP = putLogP,
                 preinitialize=preinitialize,
                 initialize=initialize,
                 # postproc=postproc,
                 tri = trind.generator(K - 1), ## symmetric indices for accessing derivative arrays
                 residuals=residuals,
                 linfo = stats,
                 rd=rd,
                 dev.resids = dev.resids,
                 linkinv = stats$linkinv, # MAYBE IT'S NEEDED IN gam.fit5
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 1, ## signal only first derivatives available...
                 discrete.ok = TRUE
  ), class = c("general.family","extended.family","family"))
} # stackFamily