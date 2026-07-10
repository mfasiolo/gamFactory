#
# A modified version of mgcv:::initial.spg as of mgcv version 1.9.4, 9th of July 2026
# Main difference marked by ### !!! ###
# Basically, for nested effects, we find lambda such that EDF is 0.5 of penalty rank.
# Instead for other effects we use mgcv:::initial.spg's cheaper heuristic.
.my.initial.spg <- function(x,y,weights,family,S,rank,off,nested_idx,offset=NULL,L=NULL,lsp0=NULL,type=1,
                            start=NULL,mustart=NULL,etastart=NULL,E=NULL,...) {
  ## initial smoothing parameter values based on approximate matching 
  ## of Frob norm of XWX and S. If L is non null then it is assumed
  ## that the sps multiplying S elements are given by L%*%sp+lsp0 and 
  ## an appropriate regression step is used to find `sp' itself.
  ## This routine evaluates initial guesses at W.
  ## Get the initial weights...
  if (length(S)==0) return(rep(0,0))
  ## start <- etastart <- mustart <- NULL
  nobs <- nrow(x) ## ignore codetools warning - required for initialization
  if (is.null(mustart)) mukeep <- NULL else mukeep <- mustart 
  eval(family$initialize) 
  if (inherits(family,"general.family")) { ## Cox, gamlss etc...   
    lbb <- family$ll(y,x,start,weights,family,offset=offset,deriv=1)$lbb ## initial Hessian
    lambda <- rep(0,length(S))
    if (TRUE) { ## experimental
      for (i in 1:length(S)) {
        ind <- off[i]:(off[i]+ncol(S[[i]])-1)
        if (rank[i]<ncol(S[[i]])) { ## find a basis for row/col space of S[[i]] and project into that.
          suppressWarnings(cs <- chol(S[[i]],pivot=TRUE))
          piv <- attr(cs,"pivot")
          Z <- S[[i]][,piv[1:rank[i]],drop=FALSE] ## basis for the space of S[[i]]
          Z <- Z/norm(Z)
          ZHZ <- -t(Z)%*%lbb[ind,ind]%*%Z
          ZSZ <- t(Z)%*%S[[i]]%*%Z
        } else { ZHZ <- -lbb[ind,ind];ZSZ <- S[[i]] }
        ### !!! START here we differ from mgcv:::initial.spg  
        if(nested_idx[i]){
          # .initialise_lambda_nested defined below
          lambda[i] <- .initialise_lambda_nested(H = ZHZ, S = ZSZ, edf = 0.5 * ncol(ZHZ)) 
        } else {
          lambda[i] <- 0.3*norm(ZHZ,"M")/norm(ZSZ,"M")
        }
        ### !!! END of modification
      }
    } else { ## original
      ## initially work out the number of times that each coefficient is penalized
      pcount <- rep(0,ncol(lbb))
      for (i in 1:length(S)) {
        ind <- off[i]:(off[i]+ncol(S[[i]])-1)
        dlb <- -diag(lbb[ind,ind,drop=FALSE])
        indp <- rowSums(abs(S[[i]]))>max(S[[i]])*.Machine$double.eps^.75 & dlb!=0
        ind <- ind[indp] ## drop indices of unpenalized
        pcount[ind] <- pcount[ind] + 1 ## add up times penalized
      }
      ## choose lambda so that corresponding elements of lbb and S[[i]]
      ## are roughly in balance...
      for (i in 1:length(S)) {
        ind <- off[i]:(off[i]+ncol(S[[i]])-1)
        lami <- 1
        #dlb <- -diag(lbb[ind,ind])
        dlb <- abs(diag(lbb[ind,ind,drop=FALSE])) 
        dS <- diag(S[[i]])
        pc <- pcount[ind]
        ## get index of elements doing any actual penalization...
        ind <- rowSums(abs(S[[i]]))>max(S[[i]])*.Machine$double.eps^.75 & dlb!=0 ## dlb > 0
        ## drop elements that are not penalizing
        dlb <- dlb[ind]/pc[ind] ## idea is to share out between penalties
        dS <- dS[ind]
        rm <- max(length(dS)/rank[i],1) ## rough correction for rank deficiency in penalty
        #while (mean(dlb/(dlb + lami * dS * rm)) > 0.4) lami <- lami*5
        #while (mean(dlb/(dlb + lami * dS * rm )) < 0.4) lami <- lami/5 
        while (sqrt(mean(dlb/(dlb + lami * dS * rm))*mean(dlb)/mean(dlb+lami*dS*rm)) > 0.4) lami <- lami*5
        while (sqrt(mean(dlb/(dlb + lami * dS * rm))*mean(dlb)/mean(dlb+lami*dS*rm)) < 0.4) lami <- lami/5
        lambda[i] <- lami 
        ## norm(lbb[ind,ind])/norm(S[[i]])
      }
    }
  } else { ## some sort of conventional regression
    if (is.null(mukeep)) {
      if (!is.null(start)) etastart <- drop(x%*%start)
      if (!is.null(etastart)) mustart <- family$linkinv(etastart)
    } else mustart <- mukeep
    if (inherits(family,"extended.family")) {
      theta <- family$getTheta()
      ## use 'as.numeric' - 'drop' can leave result as 1D array...
      Ddo <- family$Dd(y,mustart,theta,weights)
      mu.eta2 <-family$mu.eta(family$linkfun(mustart))^2 
      w <- .5 * as.numeric(Ddo$Dmu2 * mu.eta2)
      if (any(w<0)) w <- .5 * as.numeric(Ddo$EDmu2 * mu.eta2) 
    } else w <- as.numeric(weights*family$mu.eta(family$linkfun(mustart))^2/family$variance(mustart))
    w <- sqrt(w)
    if (type==1) { ## what PI would have used
      lambda <-  initial.sp(w*x,S,off)
    } else { ## balance frobenius norms
      csX <- colSums((w*x)^2) 
      lambda <- rep(0,length(S))
      for (i in 1:length(S)) {
        ind <- off[i]:(off[i]+ncol(S[[i]])-1)
        lambda[i] <- sum(csX[ind])/sqrt(sum(S[[i]]^2))
      }
    }
  }
  if (!is.null(L)) {
    lsp <- log(lambda)
    if (is.null(lsp0)) lsp0 <- rep(0,nrow(L))
    lsp <- as.numeric(coef(lm(lsp~L-1+offset(lsp0))))
    lambda <- exp(lsp)
  }
  
  return(lambda) ## initial values
  
} ## initial.spg


# Find lambda that such that trace((H+S)^-1 %*% S) matches the desired EDF
.initialise_lambda_nested <- function(H, S, edf){
  
  obj <- function(rho, edf){
    edf_rho <- sum(diag(solve(H + exp(rho) * S, H)))
    return(edf_rho - edf)
  }
  
  b <- log(c(10^-3, 10^4))
  
  for(ii in 1:10){
    lok <- obj(b[1], edf) > 0
    uok <- obj(b[2], edf) < 0
    if ( lok && uok ){
      break
    }
    if ( !lok ){
      b[1] <- b[1] + log(1/5)
    }
    if ( !uok ){
      b[2] <- b[2] + log(5)
    }
    if(ii == 10){
      stop("Impossible to initialize smoothing parameters!")
    }
  }
  
  lambda <- exp( uniroot(obj, b, tol = 0.01, edf = edf)$root )
  return(lambda)
}