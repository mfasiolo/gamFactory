#'
#' Creating a smooth constrained effect
#' 
#' @description This function creates an object representing a smooth constrained effect. 
#' At the moment, there is one effect that is constrained to be nondecreasing
#' @name createSCEffect
#' @rdname createSCEffect
#' @export createSCEffect
#'
createSCEffect <- function(x = NULL){
  
  # Loss and derivatives wrt beta (param)
  derObj <- function(param, deriv = 0) {
    
    d <- length(x)
    bt <- c(param[1], exp(param[- 1])) # beta tilde
    
    mu <- l <- drop( sum(x * bt) )
    
    if( deriv > 0 )
    {
      
      l1 <- drop( x * c(1, bt[- 1]) )
      
      if( deriv > 1 ){
        l2 <- matrix(0, d, d)
        diag(l2[- 1, - 1]) <- l1[- 1]
        
        if( deriv > 2 ){
          l3 <- list()
          coun <- 1
          for(jj in 1:d){
            for(kk in jj:d){
              for(ll in kk:d){
                l3[[coun]] <- ifelse(jj == kk & kk == ll, l2[jj, jj], 0)
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
    return( createSCEffect(x = rnorm(d)) )
  }
  
  return( list("derObj" = derObj, "initialize" = initialize) )
  
}
