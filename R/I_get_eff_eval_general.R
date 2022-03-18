.get_eff_eval_general <- function(){
  
  .eval <- function(param, deriv = 0){
    
    na <- ncol( Xi ) + 1
    nb <- length(param) - na 
    
    # Scale, inner and outer coefficients
    a0 <- exp( param[1] )
    alpha <- param[ 2:na ]
    beta <- param[ -(1:na) ]
    
    inner <- base::eval(incall)
    
    # Build P-spline basis and its derivatives
    # We are also getting the derivatives of the inner linear predictor w.r.t.
    # the scale parameter a0 and the coefficients alpha
    store <- basis$eval(x = a0 * (inner$d0 - mean(inner$d0)), deriv = deriv)
    store$g <- a0 * (inner$d0 - mean(inner$d0)) 
    store$Xi <- Xi
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$g1 <-  cbind(store$g, a0 * t(t(inner$d1) - colMeans(inner$d1)))
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        store$g2 <- cbind(store$g1, a0 * t(t(inner$d2) - colMeans(inner$d2)))
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
          store$g3 <- cbind(store$g2, a0 * t(t(inner$d3) - colMeans(inner$d3)))
        }
      }
    }
    
    o <- base::eval(efcall)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$na <- na
    o$store <- store
    o$deriv <- deriv
    
    return( o )
    
  }
  
  return( .eval )
  
}
  