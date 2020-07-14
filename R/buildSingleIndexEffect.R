#'
#' Build single index effects
#' 
#' @name buildSingleIndexEffect
#' @rdname buildSingleIndexEffect
#' @export buildSingleIndexEffect
#'
buildSingleIndexEffect <- function(Xi, splineDes){
  
  eval <- function(param, deriv = 0){
    
    na <- ncol( Xi )
    nb <- length(param) - na 
    
    # Single index and spline coefficients
    alpha <- param[ 1:na ]
    beta <- param[ -(1:na) ]
    
    # Project covariates on single index vector 
    ax <- drop( Xi %*% alpha )
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- splineDes(x = ax, deriv = deriv)
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$Xi <- Xi
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
        }
      }
    }
    
    o <- buildSingleIndexEffect(Xi = Xi, splineDes = splineDes)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$na <- na
    o$store <- store
    
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("singleIndex"))
  
  return( out )
  
}









