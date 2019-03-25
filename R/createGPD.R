#'
#' Creating a Generalized Pareto Distribution
#' 
#' @description This function creates an object representing a GPD. 
#' @name createGPDs
#' @param y a vector of observations.
#' @rdname createGPD
#' @export createGPD
#' 
createGPD <- function(y = NULL){
  
  # Likelihood and derivatives
  derObj <- function(param, deriv = 0) {
    
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    
    m <- param[ , 1, drop = TRUE]
    s <- param[ , 2, drop = TRUE]
    x <- param[ , 3, drop = TRUE]
    
    # 1) Calculate derivative of likelihood to appropriate order
    l <- - log(s) - (1 + 1/x) *  log1p( x*(-m + y)/s )
    
    if( deriv > 0 )
    {
      l_m <- (1 + x)/(s - m*x + x*y)
      l_s <- (l_m * (y-m) - 1) / s
      l_x <- - (l_m * (y-m) + (l + log(s)) / (x + 1)) / x
      
      if( deriv > 1 ){
        l_mm <- l_m^2 * x / (1+x) 
        l_ms <- - l_mm / x
        l_mx <- l_mm * (m + s - y) / (x*(1+x))
        l_ss <- - (l_ms + 1/s^2) / x
        l_sx <- l_mx * (y - m) / s
        l_xx <- - (l_x + l_mx * (y-m) + l_x / (x+1) - (l+log(s)) / (x+1)^2 ) / x
        
        if( deriv > 2 ){
          l_mmm <- 2 * l_m * l_mm * x / (1+x)
          l_mms <- 2 * l_m * l_ms * x / (1+x)
          l_mmx <- (l_mmm * (m + s - y) + l_mm)/(x*(x+1))
          l_mss <- - l_mms / x
          l_mxs <- - l_mx / s + (y-m)/s * l_mmx
          l_sss <- - (l_mss - 2/s^3) / x
          l_mxx <- - (l_mmx * (y-m) + l_mx / (x+1) - l_m / (x+1)^2 ) / x 
          l_sxx <- (y - m)/s * l_mxx
          l_ssx <- (y - m)/s * (l_mxs - l_mx/s)
          l_xxx <- - ((2*x+3)/(x+1)*l_xx + l_mxx * (y-m) - 2 * l_x / (x+1)^2 + 2*(l+log(s))/(x+1)^3 ) / x
        }
      }
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(SUM = TRUE){
      if( SUM ) { return(sum(l)) } else { return(l) }
    }
    d1 <- function(SUM = TRUE){
      if(deriv < 1) { stop("deriv < 1") }
      out <- list(l_m = l_m, l_s = l_s, l_x = l_x) 
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d2 <- function(SUM = TRUE){
      if(deriv < 2) { stop("deriv < 2") }
      out <- list(l_mm = l_mm, l_ms = l_ms, l_mx = l_mx, l_ss = l_ss, l_sx = l_sx, l_xx = l_xx)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d3 <- function(SUM = TRUE){
      if(deriv < 3) { stop("deriv < 3") }
      out <- list(l_mmm = l_mmm, l_mms = l_mms, l_mmx = l_mmx, l_mss = l_mss, l_msx = l_mxs,
                  l_mxx = l_mxx, l_sss = l_sss, l_ssx = l_ssx, l_sxx = l_sxx, l_xxx = l_xxx)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
  }
  
  # Function to generate random variables
  rd <- function(n, param) 
  {
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    mu <-  param[ , 1, drop = TRUE]
    sig <- param[ , 2, drop = TRUE]
    xi <-  param[ , 3, drop = TRUE]
    if(any(sig < 0)){ stop("scale parameter cannot be negative") }
    return( mu + sig * (runif(n)^(-xi) - 1)/xi)
  }
  
  # Random initialization of parameters
  initialize <- function(n, param){
    return( createGPD(y = rd(n, param)) )
  }
  
  return( list("derObj" = derObj, "rd" = rd, "initialize" = initialize) )
  
}
