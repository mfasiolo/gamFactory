#'
#' Build nested adaptive exponential smooth effect
#' 
#' @name eff_nexpsm
#' @rdname eff_nexpsm
#' @export eff_nexpsm
#'
eff_nexpsm <- function(y, Xi, basis, x0 = NULL){
  
  force(y); force(Xi); force(basis); force(x0)
  
  eval <- function(param, deriv = 0){
    
    na <- ncol( Xi ) + 1
    nb <- length(param) - na 
    
    # Scale, inner and outer coefficients
    a0 <- exp( param[1] )
    alpha <- param[ 2:na ]
    beta <- param[ -(1:na) ]
    
    inner <- expsmooth(y = y, Xi = Xi, beta = alpha, x0 = x0, deriv = deriv)

    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- basis(x = a0 * inner$d0, deriv = deriv)
    store$g <- a0 * inner$d0
    store$Xi <- Xi
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$g1 <-  cbind(store$g, a0 * inner$d1)
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        store$g2 <- cbind(store$g1, a0 * inner$d2)
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
          store$g3 <- cbind(store$g2, a0 * inner$d3)
        }
      }
    }
    
    o <- eff_nexpsm(y = y, Xi = Xi, basis = basis, x0 = x0)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$na <- na
    o$store <- store
    o$deriv <- deriv
    
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("nexpsm", "nested"))
  
  return( out )
  
}









