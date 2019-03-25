#'
#' Prior on direction of single index vector 
#' 
#' @description This is specific to single index effects.
#' @name createDirPrior
#' @rdname createDirPrior
#' @export createDirPrior
#'
createDirPrior <- function(r = NULL, alpha0 = NULL, alphaRho = NULL){
  
  # Loss and derivatives wrt alpha
  derObj <- function(param, deriv = 0) {
    
    d <- length( r )
    
    a <- param
    nm <- sqrt(sum(a^2))
    u <- a / nm 
    
    l <- sum( (u - r)^2 )
    
    if( deriv > 0 )
    {
      ur <- drop(crossprod(u, r))
      l1 <- drop(2 * (u * ur  - r) / nm)
      
      if( deriv > 1 ){
        l2 <- 2 * (ur * (diag(d) - 3*tcrossprod(u)) + tcrossprod(u, r) + tcrossprod(r, u)) / nm^2
        
        if( deriv > 2 ){
          l3 <- list()
          coun <- 1
          for(jj in 1:d){
            for(kk in jj:d){
              for(ll in kk:d){
                l3[[coun]] <- -3*u[ll]*l2[jj,kk]/nm +
                  2*(r[ll]*((jj==kk)-3*u[jj]*u[kk]) +
                       3*ur*(2*u[jj]*u[kk]*u[ll] - (ll==jj)*u[kk] - (ll==kk)*u[jj]) +
                       (ll==jj)*r[kk] + (ll==kk)*r[jj])/nm^3
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
  
  # Loss and derivatives wrt smoothing parameter
  derObjRho <- function(param, deriv = 0) {
    
    d <- length( r )
    
    a <- rhoAlpha(param, alpha0)$d0
    nm <- sqrt(sum(a^2))
    u <- a / nm 
    
    l <- sum( (u - r)^2 )
    ur <- drop(crossprod(u, r))
    l1 <- drop(2 * (u * ur  - r) / nm)
    l2 <- 2 * (ur * (diag(d) - 3*tcrossprod(u)) + tcrossprod(u, r) + tcrossprod(r, u)) / nm^2
    
    if( deriv ){
      
      dArho <- rhoAlpha(param, alpha0)$d1
      uTdA <- drop( crossprod(u,dArho) )
      uuT  <- tcrossprod(u,u)
      rTdA <- drop( crossprod(r,dArho) )
      uTr  <- drop( crossprod(u,r) )
      dAuT <- tcrossprod(dArho,u)
      dArT <- tcrossprod(dArho,r)
      lH <- -3*l2*uTdA/nm + 2*(rTdA*(diag(d)-3*uuT) + 3*uTr*(2*uuT*uTdA-dAuT-t(dAuT)) + dArT + t(dArT))/nm^3
      lH <- flattenSymMat( lH )
      
    }
    
    # 2) Provide accessors to derivatives
    d0 <- function(){
      return( flattenSymMat( l2 ) )
    }
    d1 <- function(){
      if(deriv < 1) { stop("deriv < 1") }
      return( lH )
    }
    
    return( list("d0" = d0, "d1" = d1) )
    
  }
  
  initialize <- function(d, alpha0 = NULL, alphaRho = NULL, ...){
    return( createDirPrior(r = rnorm(d), alpha0 = alpha0, alphaRho = alphaRho) )
  }
  
  return( list("derObj" = derObj, "derObjRho" = derObjRho, "initialize" = initialize) )
  
}
