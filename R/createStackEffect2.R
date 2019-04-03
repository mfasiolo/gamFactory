#'
#' Create derivatives of stacked effect
#' 
#' @description This is specific to stacked effects.
#' @name createStackEffect
#' @rdname createStackEffect
#' @export createStackEffect
#'
createStackEffect <- function(x = NULL){
  
  # derivatives eta (here mu) wrt nu (here param)
  derObj <- function(param, deriv = 3) {
    
    if (is.vector(param)) param <- matrix(param, nrow = 1)
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    
    d <- ncol(x) - 1
    n <- nrow(x)
    
    a <- cbind(1, exp(param)) / (1 + rowSums(exp(param)))
    am <- a[, - 1]
    xm <- x[, - 1]
    
    mu <- l <- rowSums(x * a)
    
    if( deriv > 0 )
    {
      
      l1 <- am * (xm - mu)
      if (is.numeric(am)) am <- matrix(am, nrow = n)
      if (is.numeric(l1)) l1 <- matrix(l1, nrow = n)
      
      if( deriv > 1 ){
        l2 <- list()
        coun <- 1
        # l2 <- matrix(0, d, d)
        for(jj in 1:d){
          for(kk in jj:d){
            l2[[coun]] <- l1[, jj] * (as.numeric(jj==kk) - am[, kk]) - am[, jj] * l1[, kk]
            coun <- coun + 1
            # l2[jj, kk] <- l2[kk, jj] <- l1[jj] * (as.numeric(jj==kk) - am[kk]) - am[jj] * l1[kk]
          }
        }
        l2 <- do.call(cbind, l2)
 
        i2 <- trind.generator(d)$i2
        if( deriv > 2 ){
          l3 <- list()
          coun <- 1
          for(jj in 1:d){
            for(kk in jj:d){
              for(ll in kk:d){
                l3[[coun]] <- l2[, i2[jj, ll]] * (as.numeric(kk == jj) - am[, kk]) - 
                              l1[, jj] * am[, kk] * (as.numeric(kk == ll) - am[, ll]) - 
                              l1[, kk] * am[, jj] * (as.numeric(jj == ll) - am[, ll]) - 
                              am[, jj] * l2[, i2[kk, ll]]
                coun <- coun + 1
              }
            }
          }
          l3 <- do.call(cbind, l3)
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
      if (nrow(l2) == 1) l2 <- as.numeric(l2)
      return( l2 )
    }
    d3 <- function(){
      if(deriv < 3) { stop("deriv < 3") }
      if (nrow(l3) == 1) l3 <- as.numeric(l3)
      return( l3 )
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
  }
  
  initialize <- function(n, d, ...){
    return( createStackEffect(x = matrix(rnorm(n * d), nrow = n)) )
  }
  
  return( list("derObj" = derObj, "initialize" = initialize) )
  
}