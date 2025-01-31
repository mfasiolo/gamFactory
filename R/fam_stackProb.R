#' Family for probabilistic additive stacking
#'
#' @description Creates a family to be used in mgcv::gam for probabilistic additive stacking.
#' 
#' @param logP n x K matrix of the log-predictive densities of the experts. 
#'             n is the total number of available observations, while K is the number of experts.
#'             Hence, the i-th columns contains the log-densities corresponding to the k-th expert.
#' @param ridgePen small ridge penalty on regression coefficients, useful for numerical stability.
#'                 Do not change it unless you know what you are doing.
#' @return A family that can be fitted via mgcv::gam.
#' @importFrom matrixStats rowMaxs
#' @importFrom stats make.link
#' @importFrom mgcv fix.family.link Xbd diagXVXd trind.generator Rrank gamlss.gH
#' @export
#' @name fam_stackProb
#' @rdname fam_stackProb
#' @examples 
#' library(gamFactory)
#' 
#' #### Very basic example: fitting a mixture with known mixture members (experts)
#' # Generate some data from a weighted mixture of normals
#' # The mixture weights change with x
#' set.seed(424)
#' n <- 1e3
#' x <- seq(0, pi, length.out = n)
#' w <- sin(x)
#' 
#' ii <- rbinom(n, size = 1, prob = w) 
#' y <- ii * rnorm(n, 0, 2) + (!ii) * rnorm(n, 4, 0.5)
#' plot(x, y)
#' 
#' # Normally our two experts would be fitted to some training data
#' # but here they are known and they are two normals N(0, 2) and N(4, 0.5).
#' # We evaluate their log-densities on the stacking data
#' logDen <- cbind(dnorm(y, 0, 2, log = TRUE), 
#'                 dnorm(y, 4, 0.5, log = TRUE))
#' 
#' # Estimate the mixture weights via additive stacking
#' fit <- gam(list(y ~ s(x)), 
#'            family = fam_stackProb(logDen), 
#'            data = data.frame(x=x, y=y))
#' 
#' # Compare estimated vs true weight of first expert or mixture component
#' prW <- predict(fit, type = "response", se = TRUE)
#' plot(x, prW$fit[ , 1], type = 'l', 
#'      main = "Estimated (black) and true (red) weight", 
#'      ylab = "Weight of 1st expert")
#' lines(x, prW$fit[ , 1] + 2*prW$se.fit[ , 1], lty = 2)
#' lines(x, prW$fit[ , 1] - 2*prW$se.fit[ , 1], lty = 2 )
#' lines(x, w, col = 2) # TRUTH
#' 
#' 
#' #### More realistic example: 
#' # We fit a simple and a more complex GAM model to several training data sets
#' # of increasing size. Then we simulate stacking data sets of the same sizes
#' # and we use them to estimate the weights of the two GAM models in the 
#' # stacking ensemble (a mixture of the two GAMs). We let the stacking weights
#' # depend on the sample size, because we expect that the more complex GAM should
#' # be give more weight for large data sets (because it has less bias).
#' 
#' # Generating training and stacking sets of increasing size:
#' # 60 data sets of size 100, 125 sets of size 20, ....
#' sizes <- c(rep(100, 60), rep(125, 20), rep(150, 5), 
#'            rep(175, 5), rep(200, 5), rep(300, 5),
#'            rep(400, 5), rep(600, 5), rep(800, 5))
#' 
#' ns <- length(sizes)
#' 
#' # Simulate training data sets from a standard GAM example 
#' datListTrain <- lapply(sizes, 
#'                        function( .n ){
#'                          gamSim(1, n = .n, dist = "binary", 
#'                                 scale = 0.33, verbose = FALSE)
#'                        })
#' datTrain <- do.call("rbind", datListTrain)
#' datTrain <- as.data.frame(datTrain)
#' datTrain$sizes <- rep(sizes, sizes)
#' 
#' # Simulate stacking data sets in the same way
#' datListStack <- lapply(sizes, 
#'                        function( .n ){
#'                          gamSim(1, n = .n, dist = "binary", 
#'                                 scale = 0.33, verbose = FALSE)
#'                        })
#' datStack <- do.call("rbind", datListStack)
#' datStack <- as.data.frame(datStack)
#' datStack$sizes <- rep(sizes, sizes)
#' 
#' # Estimate simple GAM on each training data set and 
#' # evaluate log-density (log-likelihood) of predictions on stacking sets
#' m1 <- lapply(1:ns, 
#'              function(.kk){
#'                .fit <- gam(y ~ s(x0, k = 3) + s(x1, k = 3) + 
#'                                s(x2, k = 3) + s(x3, k = 3), 
#'                            family=binomial, data = datListTrain[[.kk]], 
#'                            method="REML")
#'                .y <- datListStack[[.kk]]$y
#'                .p <- predict(.fit, newdata = datListStack[[.kk]], 
#'                              type = "response")
#'                .stack <- log(.p) * .y + log1p( - .p ) * ( !.y )
#'                return( list("fit" = .fit, "stack" = .stack) )
#'              })
#' 
#' # Estimate complex GAM on each training data set and 
#' # evaluate log-density (log-likelihood) of predictions on stacking sets
#' m2 <- lapply(1:ns, 
#'              function(.kk){
#'                .fit <- gam(y ~ s(x0, k = 10) + s(x1, k = 10) +
#'                                s(x2, k = 10) + s(x3, k = 10), 
#'                            family=binomial, data = datListTrain[[.kk]], 
#'                            method="REML")
#'                .y <- datListStack[[.kk]]$y
#'                .p <- predict(.fit, newdata = datListStack[[.kk]], 
#'                              type = "response")
#'                .stack <- log(.p) * .y + log1p( - .p ) * ( !.y )
#'                return( list("fit" = .fit, "stack" = .stack) )
#'              })
#' 
#' # Build matrix of log-densities on stacking set for both sets of GAMs
#' logPStack <- cbind(do.call("c", lapply(m1, "[[", "stack")), 
#'                    do.call("c", lapply(m2, "[[", "stack")))
#' 
#' # Fit additive stacking model where the experts weights depends on the 
#' # sample sizes 
#' fitStack <- gam(list(y ~ s(log(sizes), k = 7)), 
#'                  family = fam_stackProb(logP = logPStack), 
#'                  data = datStack)
#' 
#' # As expected, the weight of the second model (the more complex GAM) 
#' # increases with the size of the training set
#' prW <- predict(fitStack, type = "response", se = TRUE)
#' plot(datStack$sizes, prW$fit[ , 2], type = 'l', 
#'      main = "Estimated (black) and true (red) weight", 
#'      ylab = "Weight of 2nd expert", ylim = c(0.2, 1.1))
#' lines(datStack$sizes, prW$fit[ , 2] + 2*prW$se.fit[ , 2], lty = 2)
#' lines(datStack$sizes, prW$fit[ , 2] - 2*prW$se.fit[ , 2], lty = 2 )
#' 
fam_stackProb <- function(logP, ridgePen = 1e-5) {
  
  if (!is.null(ridgePen) && ridgePen <= 0) stop("ridgePen must be positive")
  
  link <- "identity"
  K <- ncol(logP)
  
  ## adjust log-densities for numerical stability
  X <- logP - matrixStats::rowMaxs(logP)
  expX <- exp(X)
  P <- exp(logP) ## predictive densities (not logs)
  
  # Preparing link functions
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
  
  ### Saving extra parameters in .GlobalEnv environment
  # logP
  assign(".logP", logP, envir = environment())
  getLogP <- function() get(".logP")
  putLogP <- function(.x) assign(".logP", .x, envir = environment(sys.function()))
  # P
  assign(".P", P, envir = environment())
  getP <- function() get(".P")
  putP <- function(.x) assign(".P", .x, envir = environment(sys.function()))
  # ridgePen
  assign(".ridgePen", ridgePen, envir = environment())
  getRidgePen <- function() get(".ridgePen")
  putRidgePen <- function(.x) assign(".ridgePen", .x, envir = environment(sys.function()))
  
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
    
    logP <- family$getLogP()
    P <- family$getP()
    ridgePen <- family$getRidgePen()
    
    MAT <- t(apply(logP, 1, function(x) {
      out <- rep(- 1, ncol(logP))
      out[which.max(x)] <- 1
      out
    }))
    MAT <- MAT[ , -1, drop = FALSE]
    
    ## should E be used unscaled or not?..
    use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
    
    if (is.null(start)) {
      
      jj <- attr(x,"lpi")
      
      # Objective function is penalized likelihood, 
      # mult is the multiplier of ( logP[ , k + 1] - logP[ , 1] )
      objFun <- function(mult) {
        
        if (is.list(x)) { ## discrete case
          start <- rep(0, max(unlist(jj)))
          for (k in 1:length(jj)) { ## loop over the linear predictors
            yt1 <- mult * MAT[ , k, drop = FALSE]
            R <- suppressWarnings(
              chol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,
                        qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[k]]) + 
                     crossprod(E[,jj[[k]], drop = FALSE]) +
                     if (!is.null(ridgePen)) diag(ridgePen, ncol(E[,jj[[k]], drop = FALSE])) else 0,
                   pivot=TRUE))
            Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,
                        x$qc,x$drop,lt=x$lpid[[k]])
            piv <- attr(R,"pivot")
            rrank <- attr(R,"rank")
            startji <- rep(0,ncol(R))
            if (rrank<ncol(R)) {
              R <- R[1:rrank,1:rrank]
              piv <- piv[1:rrank]
            }
            startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
            startji[!is.finite(startji)] <- 0
            start[jj[[k]]] <- startji
            
          } 
        } else {
          jj <- attr(x,"lpi")
          start <- rep(0,ncol(x))
          for (k in 1:length(jj)) { ## loop over the linear predictors
            yt1 <- mult * MAT[ , k]
            x1 <- x[,jj[[k]],drop=FALSE]
            e1 <- E[,jj[[k]],drop=FALSE] ## square root of total penalty
            if (!is.null(ridgePen)) {
              sqrt(ridgePen)
              e1rp <- matrix(0, nrow = nrow(E[,jj[[k]],drop=FALSE]),
                             ncol = ncol(E[,jj[[k]],drop=FALSE]))
              diag(e1rp) <- sqrt(ridgePen)
              e1 <- e1 + e1rp
            }
            if (use.unscaled) {
              qrx <- qr(rbind(x1,e1))
              x1 <- rbind(x1,e1)
              startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
              startji[!is.finite(startji)] <- 0
            } else startji <- penreg(x1,e1,yt1)
            start[jj[[k]]] <- startji ## copy coefficients back into overall start coef vector
          } ## lp loop
        }
        
        # Produce matrix of linear predictors
        nu <- sapply(1:length(jj), 
                     function(ii) {
                       if (is.list(x)) {
                         Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,
                             drop=x$drop,lt=x$lpid[[ii]])
                       } else {
                         x[, jj[[ii]], drop = FALSE] %*% start[jj[[ii]]]
                       }
                     })
        
        nu1 <- cbind(0, nu)
        nuCen <- nu1 - rowMaxs(nu1)
        a <- exp(nuCen) / rowSums(exp(nuCen))
        
        logLik <- log(rowSums(a * P))
        
        # Set -Inf to lowest finite logLik value
        logLik[!is.finite(logLik) & logLik < 0] <- min( logLik[is.finite(logLik)] )
        
        return( - sum(logLik) + sum( (E%*%start) ^ 2 ) )
        
      }
      
      # Optimize, then initialize using the optimal mult
      tmpOpt <- optimize(f = objFun, interval = c(0.01, 100))
      mult <- tmpOpt$minimum
      
      if (is.list(x)) { ## discrete case
        start <- rep(0,max(unlist(jj)))
        for (k in 1:length(jj)) { ## loop over the linear predictors
          yt1 <- mult * MAT[ , k]
          R <- suppressWarnings(
            chol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,
                      qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[k]]) + 
                   crossprod(E[,jj[[k]], drop = FALSE]) +
                   if (!is.null(ridgePen)) diag(ridgePen, ncol(E[,jj[[k]], drop = FALSE])) else 0,
                 pivot=TRUE))
          Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,
                      x$qc,x$drop,lt=x$lpid[[k]])
          piv <- attr(R,"pivot")
          rrank <- attr(R,"rank")
          startji <- rep(0,ncol(R))
          if (rrank<ncol(R)) {
            R <- R[1:rrank,1:rrank]
            piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
          startji[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startji
          
        } 
      } else {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        for (k in 1:length(jj)) { ## loop over the linear predictors
          yt1 <- mult * MAT[ , k]
          x1 <- x[,jj[[k]],drop=FALSE]
          e1 <- E[,jj[[k]],drop=FALSE] ## square root of total penalty
          if (!is.null(ridgePen)) {
            sqrt(ridgePen)
            e1rp <- matrix(0, nrow = nrow(E[,jj[[k]],drop=FALSE]),
                           ncol = ncol(E[,jj[[k]],drop=FALSE]))
            diag(e1rp) <- sqrt(ridgePen)
            e1 <- e1 + e1rp
          }
          if (use.unscaled) {
            qrx <- qr(rbind(x1,e1))
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0
          } else startji <- penreg(x1,e1,yt1)
          start[jj[[k]]] <- startji ## copy coefficients back into overall start coef vector
        } ## lp loop
      }
    }
  }) ## initialize
  
  
  ll <- function(y,x,coef,wt,family,offset=NULL,deriv=0,
                 d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    
    discrete <- is.list(x)
    lpi <- attr(x,"lpi") ## extract linear predictor index, in gamlss.gH it's jj
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    
    # Produce matrix of linear predictors
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
    
    if (deriv > 0) { ## grad and Hess
      ## the gradient...
      ln <- w - a[, - 1, drop = FALSE]
      
      ## the Hessian...
      lnn <- list()
      # lnn_rr <- list()
      coun <- 0
      for (jj in 1:(K - 1)) for (kk in jj:(K - 1)) {
        coun <- coun + 1
        lnn[[coun]] <- 
          ln[, jj] * (as.numeric(jj == kk) - w[, kk]) - 
          a[, jj + 1] * ln[, kk]
      }
      lnn <- do.call(cbind, lnn)
      # lnn_rr <- do.call(cbind, lnn_rr)
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
      }
      lnnn <- do.call("cbind", lnnn)
    } 
    
    if (deriv) {
      ## get the gradient and Hessian...
      ret <- gamlss.gH(x,lpi,ln,lnn,tri$i2,l3=lnnn,i3=tri$i3,l4=l4,i4=tri$i4,
                       d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l <- logLik
    
    if (!is.null(ridgePen)) {
      ret$l <- ret$l - 
        .5 * length(coef) * log(2 * pi) + 
        .5 * length(coef) * log(ridgePen) - 
        .5 * ridgePen * sum(coef ^ 2)
      
      if (deriv) {
        ret$lb <- ret$lb - ridgePen * coef
        ret$lbb <- ret$lbb - diag(ridgePen, nrow(ret$lbb))
      }
    }
    
    ret
    
  } # end ll fam_stackProb
  
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL) {
    ## optional function to give predicted values - idea is that 
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information. 
    ## if se = FALSE returns one item list containing matrix otherwise 
    ## list of two matrices "fit" and "se.fit"... 
    
    if (is.null(eta)) {
      discrete <- is.list(X)
      lpi <- attr(X,"lpi") 
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      } 
      K <- length(lpi) ## number of linear predictors
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,K)
      if (se) { 
        ve <- matrix(0,nobs,K) ## variance of eta
        ce <- matrix(0,nobs,K*(K-1)/2) ## covariance of eta_i eta_j
      } 
      for (i in 1:K) {
        if (discrete) {
          eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]]) 
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        }
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
        if (se) { ## variance and covariances for kth l.p.
          
          ve[,i] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,
                                           lt=X$lpid[[i]],rt=X$lpid[[i]]) else drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[i]]])*Xi)))
          ii <- 0
          if (i<K) for (j in (i+1):K) {
            ii <- ii + 1
            ce[,ii] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,
                                              lt=X$lpid[[i]],rt=X$lpid[[j]]) else drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[j]]])*X[,lpi[[j]]])))
          }
        }
      }
    } else { 
      se <- FALSE
    }
    gamma <- cbind(1,exp(eta))
    beta <- rowSums(gamma)
    gamma <- gamma/beta ## category probabilities
    vp <- gamma*0
    if (se) { ## need to loop to find se of probabilities...
      for (j in 1:(K+1)) {
        ## get dp_j/deta_k...
        if (j==1) dp <- -gamma[,-1,drop=FALSE]/beta else { 
          dp <- -gamma[,j]*gamma[,-1,drop=FALSE]
          dp[,j-1] <- gamma[,j]*(1-gamma[,j]) 
        }
        ## now compute variance... 
        vp[,j] <- rowSums(dp^2*ve)
        ii <- 0
        for (i in 1:K) if (i<K) for (k in (i+1):K) {
          ii <- ii + 1
          vp[,j] <- vp[,j] + 2 * dp[,i]*dp[,k]*ce[,ii] 
        }
        vp[,j] <- sqrt(pmax(0,vp[,j])) ## transform to se
      }
      return(list(fit=gamma,se.fit=vp))
    } ## if se
    list(fit=gamma)
  } ## multinom predict

  jacobian <- function(eta, jj){
    alpha <- cbind(1, exp(eta)) / rowSums(cbind(1, exp(eta)))
    K <- ncol(alpha)
    # D alpha / D eta
    DaDe <- sapply(1:(K - 1), function(.kk) {
      alpha[, jj] * (as.numeric(jj == .kk + 1) - alpha[, .kk + 1])
    })
    if(nrow(alpha) == 1) { DaDe <- matrix(DaDe, nrow = 1) }
    return(DaDe)
  }
  
  jacobian <- function(eta, jj){
    alpha <- cbind(1, exp(eta)) / rowSums(cbind(1, exp(eta)))
    K <- ncol(alpha)
    # D alpha / D eta
    DaDe <- sapply(1:(K - 1), function(.kk) {
      alpha[, jj] * (as.numeric(jj == .kk + 1) - alpha[, .kk + 1])
    })
    if(nrow(alpha) == 1) { DaDe <- matrix(DaDe, nrow = 1) }
    return(DaDe)
  }
  
  #rd <- function(mu,wt,scale) {
  
  #} ## rd
  
  # dev.resids <- function(a, b, c, d) y # MAYBE IT'S NEEDED IN gam.fit5
  
  structure(list(family="stackPredictiveFamily",ll=ll,nlp=K - 1,
                 link="identity",
                 getLogP = getLogP,
                 getP = getP,
                 getRidgePen = getRidgePen,
                 putLogP = putLogP,
                 putP = putP,
                 putRidgePen = putRidgePen,
                 preinitialize=preinitialize,
                 initialize=initialize,
                 # postproc=postproc,
                 tri = trind.generator(K - 1), ## symmetric indices for accessing deriv. arrays
                 residuals=residuals,
                 linfo = stats,
                 #rd=rd,
                 # dev.resids = dev.resids,
                 linkinv = stats$linkinv, # MAYBE IT'S NEEDED IN gam.fit5
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done, 
                 predict = predict,
                 jacobian = jacobian,
                 ls=1, ## signals that ls not needed here
                 available.derivs = 1, ## signal only first derivatives available...
                 discrete.ok = TRUE
  ), class = c("general.family","extended.family","family"))
} # stackFamily