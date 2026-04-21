#'
#' Derivatives of single index effects
#' 
#' @rdname DllkDbeta.si
#' @export DllkDbeta.si
#' @export
#'
DllkDbeta.si <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv == 0 ){ return( list() ) }
  
  if( is.null(param) ){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  }
  
  if( is.null(o$param) || !identical(param, o$param) || o$deriv < deriv){
    o <- o$eval(param = param, deriv = deriv)
  }
  
  alpha <- param[ 1:o$na ]
  beta <- param[ -(1:o$na) ]
  na <- length( alpha )
  nb <- length( beta ) 
  
  Xi_raw <- o$store$Xi
  positive_si <- isTRUE(o$xt$si$positive_si)
  a0 <- o$a0
  if(is.null(a0)) a0 <- rep(0, na)
  
  # Jacobian
  if (positive_si) {
    Xi <- t(t(Xi_raw) * exp(alpha + a0))
  } else {
    Xi <- Xi_raw
  }
  
  X <- o$store$X0; X1 <- o$store$X1; X2 <- o$store$X2; X3 <- o$store$X3
  f <- o$f
  
  der1 <- der2 <- der3 <- der4 <- NULL
  
  f1 <- drop( X1 %*% beta )
  le <- llk$d1
  lg <- le * f1
  
  ll_a <- t(Xi) %*% lg
  ll_b <- t(X) %*% le
  der1 <- c(ll_a, ll_b)
  
  if( deriv >  1){
    f2 <- drop( X2 %*% beta )
    lee <- llk$d2
    lgg <- le * f2 + lee * f1^2
    leg <- lee * f1
    
    ll_aa <- t(Xi) %*% (lgg * Xi) 
    
    # diagnoal correct for second derivative of single index with positive constraint
    if (positive_si) {
      ll_aa <- ll_aa + diag(as.vector(ll_a))
    }
    
    ll_bb <- t(X) %*% (lee * X) 
    ll_ba <- t(X) %*% (leg * Xi) + t(X1) %*% (le * Xi)
    
    der2 <- rbind(cbind(ll_aa, t(ll_ba)), 
                  cbind(ll_ba, ll_bb))
    
    der3 <- list()
    if( deriv >  2  ){
      f3 <- drop( X3%*%beta )
      leee <- llk$d3
      lggg <- le * f3 + 3 * lee * f1 * f2 + leee * f1^3
      leeg <- leee * f1
      legg <- leee * f1^2 + lee * f2
      
      coun <- 1
      for(jj in 1:na){      # AXX
        XJ <- Xi[ , jj]
        for(kk in jj:na){   # AAX
          XJK <- XJ * Xi[ , kk]
          
          for(ll in kk:na){ # AAA
            val <- sum(lggg * XJK * Xi[ , ll])
            
            # correct for third derivative of single index with positive constraint
            if (positive_si) {
              if (jj == kk && kk == ll) {
                val <- val + sum(3 * lgg * XJK + lg * XJ)
              } else if (jj == kk && kk != ll) {
                val <- val + sum(lgg * XJ * Xi[, ll])
              } else if (jj != kk && kk == ll) {
                val <- val + sum(lgg * Xi[, jj] * XJK)
              }
            }
            der3[[coun]] <- val
            coun <- coun + 1
          }
          
          for(ll in 1:nb){  # AAB
            val <- sum((legg*X[,ll]+2*leg*X1[,ll]+le*X2[,ll]) * XJK)
            # correct for AAB 
            if (positive_si && jj == kk) {
              val <- val + sum((leg * X[, ll] + le * X1[, ll]) * Xi[, jj])
            }
            der3[[coun]] <- val
            coun <- coun + 1
          }
        }
        for(kk in 1:nb){    # ABB
          XJK <- XJ * X[ , kk]
          for(ll in kk:nb){
            der3[[coun]] <- sum(XJ * (leeg*X[,ll]*X[,kk] + lee*(X1[,ll]*X[,kk]+X[,ll]*X1[,kk]) ) )
            coun <- coun + 1
          }
        }
      }
      for(jj in 1:nb){      # BBB
        XJ <- leee * X[ , jj]
        for(kk in jj:nb){
          XJK <- XJ * X[ , kk]
          for(ll in kk:nb){
            der3[[coun]] <- sum( XJK * X[ , ll] )
            coun <- coun + 1
          }
        }
      }
      der3 <- do.call("c", der3)
    } 
  }   
  
  out <- list("d1" = der1, "d2" = der2, "d3" = der3)
  return( out )
}