######
# Penalty on variance of a single index vector
#
.pen_var_si <- function(o, v, deriv = 0){
  
  na <- o$na
  a <- o$param[1:na] + o$a0
  x <- o$store$Xi
  
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

# Derivatives of penalty's Hessian w.r.t. smoothing parameters
.pen_var_si_outer <- function(o, v, DaDr){
  
  na <- o$na
  a <- o$param[1:na] + o$a0
  x <- o$store$Xi
  
  m <- ncol(DaDr)
  n <- nrow(x)  
  S <- cov(x) * (n-1) / n
  Sa <- S %*% a
  
  d1H <- list()
  for(ii in 1:m) {
    SD <- S %*% DaDr[ , ii]
    d1H[[ii]] <- 8 * ( tcrossprod(Sa, SD) + tcrossprod(SD, Sa) + S * drop(crossprod(a, S %*% DaDr[ , ii])) )
  } 

  return( d1H ) 
  
}