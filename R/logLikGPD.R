#'
#' Log-likelihood of Generalized Pareto Distribution
#' 
#' @description XXX.
#' @param np XXX.
#' @name logLikGPD
#' @rdname logLikGPD
#' @export logLikGPD
#' @examples
#' library(gamFactory)
#' n <- 1000
#' pars <- c(1, 0.3)
#' obj <- logLikGPD( )$initialize(n, pars)
#' 
#' # Derivatives should match exactly
#' fdDeriv(obj = derFunWrapper(obj$derObj), 
#'         param = pars, 
#'         ord = 1:3)
#' 
#' # Should look fine
#' der <- derivCheck(np = 100, 
#'                   parSim = function(n){ cbind(runif(n, 0.1, 1e6), runif(n, -0.45, 0.45)) }, 
#'                   obj = obj,
#'                   ord = 1:3, 
#'                   trans = function(.x){
#'                     si <- sign(.x)
#'                     return( si * sqrt(abs(.x)) )
#'                   }, 
#'                   n = 1000)
#' 
#' # We see some discrepancies here
#' par(mfrow = c(2, 2))
#' for(ii in 1:3) { plot(der[[ii]][ , 1] - der[[ii]][ , 2]) }
#' 
#' # But relative error seems acceptable 
#' par(mfrow = c(2, 2))
#' for(ii in 1:3) { plot( (der[[ii]][ , 1] - der[[ii]][ , 2]) / 
#'                          abs( der[[ii]][ , 2] )) }
#'
logLikGPD <- function(y = NULL){
  
  derObj <- function(param, deriv){
    
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    
    phi <- param[ , 1, drop = TRUE]
    xi <- param[ , 2, drop = TRUE]
    
    tol <- 1e-7
    if( any(abs(xi) < tol) ){ 
      xi[xi < tol & xi > 0] <- tol
      xi[xi > -tol & xi <= 0] <- - tol
    }
    
    Cxi <- xi * (1 + xi) * y
    
    llO <- log1p(xi) - log(phi) - (1/xi+1) * log1p(Cxi / phi)
    
    if( deriv > 0 )
    {
      
      llO_xi <- 1/(1+xi)+1/xi^2*log1p(Cxi/phi) - (1/xi+1)*(2*xi+1)*y/(phi+Cxi)
      llO_phi <- (-phi+y+xi*y)/(phi*(phi+Cxi))
      
      if( deriv > 1 ){
        
        llO_xi2 <- - 1/(1+xi)^2-2/xi^3*log1p(Cxi/phi) + 
          2/xi^2*(2*xi+1)*y/(phi+Cxi) -
          (1/xi+1)*(2*y*(phi+Cxi)-(2*xi+1)^2*y^2)/(phi+Cxi)^2
        
        llO_phi2<- (phi^2-2*phi*(1+xi)*y-xi*(1+xi)^2*y^2)/(phi^2*(phi+Cxi)^2)
        
        llO_phixi<- -(1+xi)*y*(-2*phi+y+xi*y)/(phi*(phi+Cxi)^2)
        
        if( deriv > 2 ){
          
          A <- phi+xi*(1+xi)*y
          term_1 <- 2/(1+xi)^3-2*(1+xi)*(y+2*xi*y)^3/(xi*A^3)
          term_2 <- 6*(1+xi)*(1+2*xi)*y^2/(xi*A^2)-3*(y+2*xi*y)^2/(xi^2*A^2)
          term_3 <- 6*y/(xi^2*A)-6*(y+2*xi*y)/(xi^3*A)
          term_4 <- 6*log1p(Cxi/phi)/xi^4
          llO_xi3<- term_1 + term_2 + term_3 + term_4
          
          llO_phi3<- (-2*phi^3+6*phi^2*(1+xi)*y+6*phi*Cxi^2/xi+2*Cxi^3/xi)/(phi^3*(phi+Cxi)^3)
          
          llO_phixi2<- 2*y*(phi^2-3*phi*(1+xi)^2*y+(1+xi)^3*y^2)/(phi*(phi+Cxi)^3)
          
          llO_phi2xi<-(1+xi)*y*(-4*phi^2+3*phi*(1+xi)*y+xi*(1+xi)^2*y^2)/(phi^2*(phi+Cxi)^3)
          
        }
      }
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(SUM = TRUE){
      if( SUM ) { return(sum(llO)) } else { return(llO) }
    }
    d1 <- function(SUM = TRUE){
      if(deriv < 1) { stop("deriv < 1") }
      out <- list(llO_phi=llO_phi, llO_xi = llO_xi) 
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d2 <- function(SUM = TRUE){
      if(deriv < 2) { stop("deriv < 2") }
      out <- list(llO_phi2=llO_phi2, llO_phixi=llO_phixi, llO_xi2=llO_xi2)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d3 <- function(SUM = TRUE){
      if(deriv < 3) { stop("deriv < 3") }
      out <- list(llO_phi3=llO_phi3, llO_phi2xi=llO_phi2xi, llO_phixi2=llO_phixi2, llO_xi3=llO_xi3)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
    
  }
  
  rd <- function(n, param) 
  {
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    phi <-  param[1]
    xi <- param[2]
    sig <- phi / (1 + xi)
    if(any(sig < 0)){ stop("scale parameter cannot be negative") }
    return( sig * ((1-runif(n))^(-xi) - 1)/xi)
  }
  
  initialize<-function(n, param){
    return( logLikGPD(y = rd(n, param)) )
  }
  
  return( list("derObj" = derObj, "rd" = rd, "initialize" = initialize) )
  
  
}
