#'
#' Log-likelihood of a Markov process
#' 
#' @description XXX.
#' @param np XXX.
#' @name logLikMARKOV
#' @rdname logLikMARKOV
#' @export logLikMARKOV
#' @examples
#' 
#' library(gamFactory)
#' 
#' # Simulate some data from Markov chain with two states (0, 1)
#' n <- 5000
#' x = sort(runif(n, -1, 1))
#' 
#' # Probabilities of staying in state 0 and 1
#' p1 <- exp(x^2) / (1 + exp(x^2))
#' p2 <- exp(-x - 0.6*x^2) / (1 + exp(-x - 0.6*x^2))
#' 
#' # Simulate observations 
#' y0 <- numeric(n)
#' y0[1] <- 1
#' for(ii in 2:n){
#'   if(y0[ii-1] == 0) y0[ii] <- rbinom(1, 1, 1 - p1[ii])
#'   if(y0[ii-1] == 1) y0[ii] <- rbinom(1, 1, p2[ii])
#' }
#' 
#' myDat <- data.frame(x = x, y = y0)
#' 
#' # Plot data and probabilities
#' par(mfrow = c(2, 2))
#' plot(myDat)
#' plot(x, p1)
#' plot(x, p2)
#' 
#' # Create Markov model and fit it
#' library(mgcViz)
#' myMRK <- buildGAMLSS( bundleMARKOV )
#' fit <- gamV(list(y ~ s(x), ~ s(x)), data = myDat, family = myMRK)
#' 
#' # Compare true and estimated probabilities
#' plot(fit, select = 1, 
#'      trans = function(.x) exp(coef(fit)[1] + .x) / (1 + exp(coef(fit)[1] + .x))) + 
#'   l_fitLine() + l_ciLine() + 
#'   geom_line(data = data.frame(x = x, y = p1), 
#'             mapping = aes(x=x, y=y), inherit.aes = FALSE, color = "red")
#' 
#' plot(fit, select = 2, 
#'      trans = function(.x) exp(coef(fit)[11] + .x) / (1 + exp(coef(fit)[11] + .x))) + 
#'   l_fitLine() + l_ciLine() + 
#'   geom_line(data = data.frame(x = x, y = p2), 
#'             mapping = aes(x=x, y=y), inherit.aes = FALSE, color = "red")
#'
logLikMARKOV <- function(y = NULL){
  
  derObj <- function(param, deriv){
    
    ly <- y[ -length(y) ]
    dy <- diff( y )
    n <- length( dy )
    t11 <- which( ly == 0 & dy == 0 )
    t22 <- which( ly == 1 & dy == 0 )
    t12 <- which( dy == 1 )
    t21 <- which( dy == -1 )
    
    if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
    
    p1 <- param[ , 1, drop = TRUE]
    p2 <- param[ , 2, drop = TRUE]
    
    if( length(p1) == 1 ){
      p1 <- rep(p1, n)
      p2 <- rep(p2, n)
    } else {
      p1 <- p1[ -1 ]
      p2 <- p2[ -1 ]
    }
  
    ll_0 <- numeric( n )
    ll_0[t11] <- log(p1[t11])
    ll_0[t22] <- log(p2[t22])
    ll_0[t12] <- log1p(-p1[t12])
    ll_0[t21] <- log1p(-p2[t21])
    ll_0 <- c(0, ll_0)
    
    if( deriv > 0 )
    {
      ll_1 <- ll_2 <- numeric( n )
      ll_1[t11] <- 1 / p1[t11]
      ll_2[t22] <- 1 / p2[t22]
      ll_1[t12] <- - 1 / (1 - p1[t12])
      ll_2[t21] <- - 1 / (1 - p2[t21])
      ll_1 <- c(0, ll_1)
      ll_2 <- c(0, ll_2)
      
      if( deriv > 1 ){
        
        ll_11 <- ll_22 <- ll_12 <- numeric( n )
        ll_11[t11] <- - 1 / p1[t11]^2
        ll_22[t22] <- - 1 / p2[t22]^2
        ll_11[t12] <- - 1 / (1 - p1[t12])^2
        ll_22[t21] <- - 1 / (1 - p2[t21])^2
        ll_11 <- c(0, ll_11)
        ll_22 <- c(0, ll_22)
        ll_12 <- c(0, ll_12)
        
        if( deriv > 2 ){
          
          ll_111 <- ll_112 <- ll_122 <- ll_222 <- numeric( n )
          ll_111[t11] <- 2 / p1[t11]^3
          ll_222[t22] <- 2 / p2[t22]^3
          ll_111[t12] <- - 2 / (1 - p1[t12])^3
          ll_222[t21] <- - 2 / (1 - p2[t21])^3
          ll_111 <- c(0, ll_111)
          ll_222 <- c(0, ll_222)
          ll_112 <- c(0, ll_112)
          ll_122 <- c(0, ll_122)
          
        }
      }
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(SUM = TRUE){
      if( SUM ) { return(sum(ll_0)) } else { return(ll_0) }
    }
    d1 <- function(SUM = TRUE){
      if(deriv < 1) { stop("deriv < 1") }
      out <- list(ll_1=ll_1, ll_2=ll_2) 
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d2 <- function(SUM = TRUE){
      if(deriv < 2) { stop("deriv < 2") }
      out <- list(ll_11=ll_11, ll_12=ll_12, ll_22=ll_22)
      if( SUM ){
        out <- sapply(out, sum)
      } 
      return(out)
    }
    d3 <- function(SUM = TRUE){
      if(deriv < 3) { stop("deriv < 3") }
      out <- list(ll_111=ll_111, ll_112=ll_112, ll_122=ll_122, ll_222=ll_222)
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
    p1 <- param[ , 1]
    p2 <- param[ , 2]
    out <- numeric(n)
    out[1] <- 1
    for(ii in 2:n){
      if(out[ii-1] == 0) out[ii] <- rbinom(1, 1, 1 - p1)
      if(out[ii-1] == 1) out[ii] <- rbinom(1, 1, p2)
    }
    return( out )
  }
  
  initialize <- function(n, param){
    return( logLikMARKOV(y = rd(n, param)) )
  }
  
  return( list("derObj" = derObj, "rd" = rd, "initialize" = initialize) )
  
  
}
