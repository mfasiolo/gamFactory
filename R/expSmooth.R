#'
#' Exponential smooth and its derivatives
#' 
#' @name expSmooth
#' @rdname expSmooth
#' @export expSmooth
#'
expSmooth <- function(x, alpha, x0 = NULL, deriv = 0){
  
  if( is.matrix(x) ){ x <- as.vector(x) }
    
  if( is.null(x0) ){ x0 <- x[1] }
  
  n <- length(x)
  
  xh <- numeric(n)
  xh[1] <- alpha * x0 + (1-alpha) * x[1]
  for(ii in 2:n){
    xh[ii] <- alpha * xh[ii-1] + (1-alpha) * x[ii]
  }
  
  d1 <- d2 <- d3 <- NULL
  if( deriv ){
    d1 <- numeric(n)
    d1[1] <- -x[1] + x0
    for(ii in 2:n){
      d1[ii] <- -x[ii] + xh[ii-1] + alpha*d1[ii-1]
    }
    
    if( deriv > 1 ){
      d2 <- numeric(n)
      d2[1] <- 0
      for(ii in 2:n){
        d2[ii] <- 2*d1[ii-1] + alpha*d2[ii-1]
      }
      
      if( deriv > 2 ){
        d3 <- numeric(n)
        d3[1] <- 0
        for(ii in 2:n){
          d3[ii] <- 3*d2[ii-1] + alpha*d3[ii-1]
        }
      }
    }
  }
  
  return( list("d0" = xh, "d1" = d1, "d2" = d2, "d3" = d3) )
    
}
 