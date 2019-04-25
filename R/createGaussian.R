#'
#' Creating a Gaussian distribution
#' 
#' @description This function creates an object representing a Gaussian distribution 
#' @name createGaussian
#' @param y a vector of observations.
#' @rdname createGaussian
#' @export createGaussian
#' 
createGaussian <- function(y = NULL) {
  
  derObj <- function(param, deriv = 0) {
    
    if (is.vector(param)) param <- matrix(param, nrow = 1)
    if (ncol(param) != 2) stop("Wrong number of parameters provided")
    
    mu <- param[ , 1, drop = TRUE]
    g <- param[ , 2, drop = TRUE] ## g is log(sigma)
    tau2 <- exp( - 2 * g )
    n <- length(y)

    if (length(mu) == 1) {
      mu <- rep(mu, n)
      g <- rep(g, n)
      tau2 <- rep(tau2, n)
    }
    
    ymu <- y - mu          ## precompute
    ymu2 <- ymu ^ 2        ## precompute
    
    l <- - .5 * log(2 * pi) - g - .5 * tau2 * ymu2
    
    if( deriv > 0 )
    {

      ## First derivatives 
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
    if (ncol(param != 2)) stop("Wrong number of parameters provided")
    mu <-  param[ , 1, drop = TRUE]
    sig <- exp( param[ , 2, drop = TRUE] )
    
    return( rnorm(n, mean = mu, sd = sig) )
  }
  
  initialize <- function(n, param, ...) {
    return( createGaussian(y = rd(n = n, param = param)) )
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

