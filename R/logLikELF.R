#'
#' Log-likelihood of Extended Log-F (ELF) distribution
#' @description XXX.
#' @param np XXX.
#' @name logLikELF
#' @rdname logLikELF
#' @export logLikELF
#' @examples 
#' library(gamFactory)
#' n <- 1000
#' pars <- c(1, 0.5)
#' obj <- logLikELF( )$initialize(n, pars, h = 0.7, qu = 0.4)
#' 
#' # Derivatives should match exactly
#' fdDeriv(obj = derFunWrapper(obj$derObj), 
#'         param = pars, 
#'         ord = 1:3)
#' 
#' # Should look fine
#' der <- derivCheck(np = 100, 
#'                   parSim = function(n){ cbind(rnorm(n, 0, 1e3), 0.1 + rexp(n, 1)) }, 
#'                   obj = obj,
#'                   ord = 1:3, 
#'                   trans = function(.x){
#'                     si <- sign(.x)
#'                     return( si * sqrt(abs(.x)) )
#'                   }, 
#'                   qu = 0.99, 
#'                   h = 0.2, 
#'                   n = 1000)
#' 
#' par(mfrow = c(2, 2))
#' for(ii in 1:3) { plot(der[[ii]][ , 1] - der[[ii]][ , 2]) }
#' 
#' par(mfrow = c(2, 2))
#' for(ii in 1:3) { plot( (der[[ii]][ , 1] - der[[ii]][ , 2]) / 
#'                          abs( der[[ii]][ , 2] )) }
#' 
logLikELF <- function(y, qu, h){
  
  derObj <- function(param, deriv){
    
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    
    n <- nrow( param ) 
    
    tau <- qu
    mu <- param[ , 1, drop = TRUE]
    sig <- param[ , 2, drop = TRUE]
    lam <- h / sig
    
    z <- (y - mu) / sig
    zc <- (y - mu) / h
    lpxp <- log1pexp( zc ) 
    
    ll0 <- (1-tau) * z - h*lpxp/sig - log( h*beta(h*(1-tau)/sig, h*tau/sig) )
    
    if( deriv > 0 )
    {
      # D logBeta / D lam;  D^2 logBeta / D lam^2 
      dBdL <- (1-tau) * digamma(lam*(1-tau)) + tau * digamma(lam*tau) - digamma(lam)  
      d2BdL2 <- (1-tau)^2 * trigamma(lam*(1-tau)) + tau^2 * trigamma(lam*tau) - trigamma(lam)  
      
      # D lam / D sig;  D^2 lam / D sig^2 
      dLdS <- - h / sig^2
      d2LdS2 <- 2 * h / sig^3
      
      # D logBeta / D sig;  D^2 logBeta / D sig^2
      dBdS <- dLdS * dBdL
      d2BdS2 <- d2LdS2*dBdL + (dLdS)^2*d2BdL2
      
      gLog <- (tau-1)*(y-mu) + h * lpxp
      
      dl <- dlogis(y, mu, h)
      pl <- plogis(y, mu, h)
      
      dLLKdmu <- (pl - 1 + tau) / sig
      ll_m <- dLLKdmu
      ll_s <- gLog/sig^2 - dBdS
   
      if( deriv > 1 ){
        
        D2LLKdmu2 <- - dl / sig
        l2 <- matrix(0, n, 3)
        ## order mm,ms,ss
        ll_mm <- D2LLKdmu2
        ll_ms <- - dLLKdmu / sig
        ll_ss <- - 2*gLog/sig^3 - d2BdS2
        
        if( deriv > 2 ){
          
          # D^3 logBeta / D lam^3 ;  D^3 lam / D sig^3;  D^3 logBeta / D sig^3
          d3BdL3 <- (1-tau)^3 * psigamma(lam*(1-tau), 2) + tau^3 * psigamma(lam*tau, 2) - psigamma(lam, 2)
          d3LdS3 <- - 6 * h / sig^4
          d3BdS3 <- d3LdS3*dBdL + 3*dLdS*d2LdS2*d2BdL2 + dLdS^3*d3BdL3
          
          der <- sigmoid(zc, deriv = TRUE)
          
          D3LLKdmu3 <- der$D2 / (sig*h^2)
          
          ## the third derivatives
          ## order mmm,mms,mss,sss
          ll_mmm <- D3LLKdmu3
          ll_mms <- - D2LLKdmu2 / sig
          ll_mss <- 2 * dLLKdmu / sig^2
          ll_sss <- 6*gLog/sig^4 - d3BdS3
   
          
          if( deriv > 3 ){
            # D^4 logBeta / D lam^4 ;  D^4 lam / D sig^4;  D^4 logBeta / D sig^4
            d4BdL4 <- (1-tau)^4 * psigamma(lam*(1-tau), 3) + tau^4 * psigamma(lam*tau, 3) - psigamma(lam, 3)
            d4LdS4 <- 24 * h / sig^5
            d4BdS4 <- d4LdS4*dBdL + 3*d2LdS2^2*d2BdL2 + 4*dLdS*d3LdS3*d2BdL2 + 6*(dLdS)^2*d2LdS2*d3BdL3 + dLdS^4*d4BdL4 
            
            ## the fourth derivatives, order: mmmm,mmms,mmss,msss,ssss
            ll_mmmm <- - der$D3 / (sig*h^3)
            ll_mmms <- - D3LLKdmu3 / sig 
            ll_mmss <- 2 * D2LLKdmu2 / sig^2
            ll_msss <- - 6 * dLLKdmu / sig^3
            ll_ssss <- - 24*gLog/sig^5 - d4BdS4
            
          }
          
        }
      }
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(SUM = TRUE){
      if( SUM ) { return(sum(ll0)) } else { return(ll0) }
    }
    d1 <- function(SUM = TRUE){
      if(deriv < 1) { stop("deriv < 1") }
      out <- list(ll_m = ll_m, ll_s = ll_s) 
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d2 <- function(SUM = TRUE){
      if(deriv < 2) { stop("deriv < 2") }
      out <- list(ll_mm = ll_mm, ll_ms = ll_ms, ll_ss = ll_ss)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d3 <- function(SUM = TRUE){
      if(deriv < 3) { stop("deriv < 3") }
      out <- list(ll_mmm = ll_mmm, ll_mms = ll_mms, ll_mss = ll_mss, ll_sss = ll_sss)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
    
  }
  
  rd <- function(n, param, ...) 
  { # Simulate from Gaussian with random mean and variance
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    mu <-  param[1]
    sig <- rexp(1, 1)
    if(any(sig < 0)){ stop("scale parameter cannot be negative") }
    return( rnorm(n, mu, sig) )
  }
  
  initialize<-function(n, param, qu, h){
    return( logLikELF(y = rd(n, param), qu, h) )
  }
  
  ml <- function(y) {
    # For inizialization of theta parameters
    muHat <- quantile(y, qu)
    sigHat <- sd(y)
    return(c(muHat, sigHat))
  }
  
  return( list("derObj" = derObj, "rd" = rd, "initialize" = initialize, 
               "ml" = ml, npar = 2) )
  
  
}
