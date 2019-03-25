#'
#' Create derivatives of stacked effect
#' 
#' @description This is specific to stacked effects.
#' @name createStackEffect
#' @rdname createStackEffect
#' @export createStackEffect
#'
createStackEffect <- function(x = NULL){
  
  # Loss and derivatives wrt alpha
  derObj <- function(param, deriv = 0) {
    
    d <- length(x) - 1
    
    a <- c(1, exp( param )) / (1 + sum(exp(param)))
    am <- a[-1]
    xm <- x[-1]
    
    mu <- l <- drop( sum(x * a) )
    
    if( deriv > 0 )
    {
      
      l1 <- drop( am * (xm - mu) )
      
      if( deriv > 1 ){
        l2 <- matrix(0, d, d)
        for(jj in 1:d){
          for(kk in jj:d){
            l2[jj, kk] <- l2[kk, jj] <- l1[jj] * (as.numeric(jj==kk) - am[kk]) - am[jj] * l1[kk]
          }
        }
        
        if( deriv > 2 ){
          l3 <- list()
          coun <- 1
          for(jj in 1:d){
            for(kk in jj:d){
              for(ll in kk:d){
                l3[[coun]] <- l2[jj,ll] * (as.numeric(kk==jj) - am[kk]) - 
                              l1[jj] * am[kk] * (as.numeric(kk==ll) - am[ll]) - 
                              l1[kk] * am[jj] * (as.numeric(jj==ll) - am[ll]) - am[jj] * l2[kk,ll]
                coun <- coun + 1
              }
            }
          }
          l3 <- do.call("c", l3)
        }
      }
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(){
      return( l )
    }
    d1 <- function(){
      if(deriv < 1) { stop("deriv < 1") }
      return( l1 )
    }
    d2 <- function(){
      if(deriv < 2) { stop("deriv < 2") }
      return( l2 )
    }
    d3 <- function(){
      if(deriv < 3) { stop("deriv < 3") }
      return( l3 )
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
  }
  
  initialize <- function(d, ...){
    return( createStackEffect(x = rnorm(d)) )
  }
  
  return( list("derObj" = derObj, "initialize" = initialize) )
  
}