#'
#' Log-likelihood of the Gaussian distribution
#' @description XXX.
#' @param np XXX.
#' @name logLikGAU
#' @rdname logLikGAU
#' @export logLikGAU
#' @examples 
#' library(gamFactory)
#' n <- 1000
#' pars <- c(1, 0.5) # mu, log(sigma)
#' obj <- logLikGAU( )$initialize(n, pars)
#' 
#' # Derivatives should match exactly
#' fdDeriv(obj = derFunWrapper(obj$derObj), 
#'         param = pars, 
#'         ord = 1:3)
#' 
#' # Should look fine
#' der <- derivCheck(np = 100, 
#'                   parSim = function(n){ cbind(rnorm(n, 0, 1e3), log(1 + rexp(n, 1))) }, 
#'                   obj = obj,
#'                   ord = 1:3, 
#'                   trans = function(.x){
#'                     si <- sign(.x)
#'                     return( si * sqrt(abs(.x)) )
#'                   }, 
#'                   n = 1000)
#' 
#' par(mfrow = c(2, 2))
#' for(ii in 1:3) { plot(der[[ii]][ , 1] - der[[ii]][ , 2]) }
#' 
#' par(mfrow = c(2, 2))
#' for(ii in 1:3) { plot( (der[[ii]][ , 1] - der[[ii]][ , 2]) / 
#'                          abs( der[[ii]][ , 2] )) }
#' 
#' 
logLikGAU <- function(y = NULL) {
  
  derObj <- function(param, deriv = 0) {
    
    if (is.vector(param)) param <- matrix(param, nrow = 1)
    if (ncol(param) != 2) stop("Wrong number of parameters provided")
    
    mu <- param[ , 1, drop = TRUE]
    g <- param[ , 2, drop = TRUE] # log(sigma)
    tau2 <- exp( - 2 * g ) # 1 / sigma^2
    n <- length(y)
    
    if (length(mu) == 1) {
      mu <- rep(mu, n)
      g <- rep(g, n)
      tau2 <- rep(tau2, n)
    }
    
    ymu <- y - mu        
    ymu2 <- ymu ^ 2   
    
    l <- - .5 * log(2 * pi) - g - .5 * tau2 * ymu2
    
    if( deriv > 0 )
    {
      
      ## First derivatives wrt mu and log(sigma)
      Dm <- tau2 * ymu
      Ds <- - 1 + tau2 * ymu2
      
      if( deriv > 1 ){
        Dmm <- - tau2
        Dms <- - 2 * tau2 * ymu
        Dss <- - 2 * tau2 * ymu2
        
        if( deriv > 2 ){
          Dmmm <- rep(0, n)
          Dmms <- 2 * tau2
          Dmss <- 4 * tau2 * ymu
          Dsss <- 4 * tau2 * ymu2
        }
      }
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(SUM = TRUE){
      if( SUM ) { return(sum(l)) } else { return(l) }
    }
    d1 <- function(SUM = TRUE){
      if(deriv < 1) { stop("deriv < 1") }
      out <- list(l_m = Dm, l_s = Ds) 
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d2 <- function(SUM = TRUE){
      if(deriv < 2) { stop("deriv < 2") }
      out <- list(l_mm = Dmm, l_ms = Dms, l_ss = Dss)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d3 <- function(SUM = TRUE){
      if(deriv < 3) { stop("deriv < 3") }
      out <- list(l_mmm = Dmmm, l_mms = Dmms, l_mss = Dmss, l_sss = Dsss)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
  }
  
  rd <- function(n, param) 
  {
    if (is.vector(param)) param <- matrix(param, nrow = 1)
    if (ncol(param) != 2) stop("Wrong number of parameters provided")
    mu <-  param[ , 1, drop = TRUE]
    sig <- exp( param[ , 2, drop = TRUE] )
    
    return( rnorm(n, mean = mu, sd = sig) )
  }
  
  initialize <- function(n, param, ...) {
    return( logLikGAU(y = rd(n = n, param = param)) )
  }
  
  ml <- function(y) {
    n <- length(y)
    muHat <- mean(y)
    logSigmaHat <- log(sd(y) * (n - 1) / n)
    return(c(muHat, logSigmaHat))
  }
  
  return( list("derObj" = derObj, "rd" = rd, "initialize" = initialize,
               ml = ml, npar = 2) )
  
}

