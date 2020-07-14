#'
#' Penalty on variance of a single index vector
#' 
#' @description Penalty on variance of a single index vector. 
#' @param a single index parameter vector of dimension p.
#' @param x an n by p matrix of covariates to be projected.
#' @param v target variance of the projected covariates.
#' @param deriv number of derivatives to be calculated.
#' @param DaDr derivative for single index vector w.r.t. smoothing parameter.
#' @name pen_varSI
#' @rdname pen_varSI
#' @export pen_varSI
#'
pen_varSI <- function(a, x, v, deriv = 0){
  
  p <- length(a)
  ax <- x %*% a
  
  # Empirical variance
  vhat <- sum(ax^2)/length(ax) - mean(ax)^2
  
  # Loss
  l0 <- (vhat - v)^2 
  
  l1 <- l2 <- l3 <- NULL
  if( deriv ){
    n <- nrow(x)
    
    S <- cov(x) * (n-1) / n
    Sa <- S %*% a
    
    l1 <- 4 * (vhat - v) * Sa
    
    if( deriv > 1 ){
      l2 <- 8 * tcrossprod(Sa, Sa) + 4 * (vhat - v) * S 
    }
    
    if( deriv > 2 ){
      for(i1 in 1:p){
        for(i2 in i1:p){
          for(i3 in i2:p){
            tmp <- 8 * (S[i1, i3] * Sa[i2] + Sa[i1] * S[i2, i3] + S[i1, i2] * Sa[i3])
            l3 <- c(l3, tmp)
          }
        }
      }
      
    }
  } 
  
  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = l3) )
  
}


#' @name pen_varSI
#' @rdname pen_varSI
#' @export pen_varSI_outer
pen_varSI_outer <- function(a, x, DaDr){
  
  n <- nrow(x)  
  S <- cov(x) * (n-1) / n
  Sa <- S %*% a
  SD <- S %*% DaDr
  
  out <- 8 * { tcrossprod(Sa, SD) + tcrossprod(SD, Sa) + S * drop(crossprod(a, S %*% DaDr)) }
  
  return( out ) 
  
}