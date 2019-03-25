#'
#' Prior on variance of single index model 
#' 
#' @description This is specific to single index effects.
#' @name createVarPrior
#' @rdname createVarPrior
#' @export createVarPrior
#'
createVarPrior <- function(X = NULL, v = NULL, alpha0 = NULL, rhoAlpha = NULL){
  
  # Derivative of penalty wrt alpha
  derObj <- function(param, deriv = 0) {
    
    d <- ncol(X)
    n <- nrow(X)
    
    al <- param
    xa <- X %*% al
    vr <- mean(xa^2) - mean(xa)^2 # Variance (denominator n, not n-1)
    
    l <- (vr - v)^2
    
    if( deriv > 0 )
    {
      S <- cov(X) * (n-1) / n
      l1 <- drop( 4 * (vr - v) * (S %*% al) )
      
      if( deriv > 1 ){
        l2 <- 8 * tcrossprod(drop(S %*% al)) + 4 * (vr - v) * S 
        
        if( deriv > 2 ){
          l3 <- list()
          coun <- 1
          for(jj in 1:d){
            for(kk in jj:d){
              for(ll in kk:d){
                l3[[coun]] <- 8 * (S[jj,ll] * crossprod(al, S[ , kk]) + 
                                     S[kk,ll] * crossprod(al, S[ , jj])) + 8 * S[jj,kk] * crossprod(al, S[ , ll]) 
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
  
  # Derivative of penalty wrt alpha
  derObjRho <- function(param, deriv = 0) {
    
    param <- drop( param )
    
    d <- ncol(X)
    n <- nrow(X)
    
    al <- rhoAlpha(param, alpha0)$d0
    xa <- X %*% al
    vr <- mean(xa^2) - mean(xa)^2 # Variance (denominator n, not n-1)
    
    l <- (vr - v)^2
    
    S <- cov(X) * (n-1) / n
    
    l2 <- flattenSymMat(8 * tcrossprod(drop(S %*% al)) + 4 * (vr - v) * S)
    
    if( deriv > 0 )
    {
      
      dArho <- rhoAlpha(param, alpha0)$d1
      lH <- 8 * (S%*%al %*% t(S%*%dArho) + S%*%dArho %*% t(S%*%al) + S*drop((t(al)%*%S)%*%dArho))
      lH <- flattenSymMat( lH )
      
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(){
      return( l2 )
    }
    d1 <- function(){
      if(deriv < 1) { stop("deriv < 1") }
      return( lH )
    }
    
    return( list("d0" = d0, "d1" = d1) )
    
  }
  
  initialize <- function(X, v, alpha0 = NULL, rhoAlpha = NULL, ...){
    return( createVarPrior(X = X, v = v, alpha0 = alpha0, rhoAlpha = rhoAlpha) )
  }
  
  return( list("derObj" = derObj, "derObjRho" = derObjRho, "initialize" = initialize) )
  
}

