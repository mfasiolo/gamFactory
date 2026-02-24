.get_eff_eval_si_nexp <- function(){
  
  .eval <- function(param, deriv = 0){
    
    # param = [alpha_nexp, alpha_si, beta]
    # prepare for inner <- base::eval(incall)
    alpha_nexp <- param[1:n_nexp]
    alpha_si <- param[(n_nexp+1):na]
    beta <- param[(na+1):length(param)]
    
    inner <- eval(incall)
    
    
    # Build P-spline basis and its derivatives
    # We are also getting the derivatives of the inner linear predictor w.r.t.
    # the scale parameter alpha_scale and the coefficients alpha_si, alpha_nexp
    store <- basis$evalX(x = (inner$d0 - mean(inner$d0)), deriv = deriv)
    store$g <- inner$d0 - mean(inner$d0)
    store$X_nexp <- X_nexp
    store$X_si <- X_si
    
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$g1 <-  t(t(inner$d1) - colMeans(inner$d1))
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        store$g2 <- t(t(inner$d2) - colMeans(inner$d2))
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
          store$g3 <- t(t(inner$d3) - colMeans(inner$d3))
        }
      }
    }
    
    o <- base::eval(efcall)
    o$f <- drop( store$X0 %*% beta )
    
    # param = [alpha_si,alpha_nexp(with alpha_scale),beta]
    o$param <- param
    o$na <- na         #length of alpha_si + alpha_nexp (without alpha_scale)
    o$store <- store
    o$deriv <- deriv
    o$alpha_center <- alpha_center
    o$positive_si <- positive_si
    
    return( o )
    
  }
  
  return( .eval )
  
}
