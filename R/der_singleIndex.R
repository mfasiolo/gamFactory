#'
#' Derivatives of single index effects
#' 
#' @rdname der.singleIndex
#' @export der.singleIndex
#' @export
#'
der.singleIndex <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv == 0 ){ return( list() ) }
  
  if( is.null(param) ){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  }
  
  # Need to update the object
  if( is.null(o$param) || !identical(param, o$param)){
    o <- o$eval(param = param, deriv = deriv)
  }
  
  # Extract stuff from object for convenience
  alpha <- param[ 1:o$na ]
  beta <- param[ -(1:o$na) ]
  
  na <- length( alpha )
  nb <- length( beta ) 
  
  Xi <- o$store$Xi
  
  # Outer design matrix and its derivatives
  X <- o$store$X0; X1 <- o$store$X1; X2 <- o$store$X2; X3 <- o$store$X3
  
  f <- o$f
  
  # Project covariates on single index vector 
  ax <- drop( Xi %*% alpha )
  
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
    ll_bb <- t(X) %*% (lee * X)    # Same as t(X) %*% diag(le2) %*% X
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
            der3[[coun]] <- sum(lggg * XJK * Xi[ , ll])
            coun <- coun + 1
          }
          for(ll in 1:nb){  # AAB
            der3[[coun]] <- sum((legg*X[,ll]+2*leg*X1[,ll]+le*X2[,ll]) * XJK)
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
            der3[[coun]] <- sum( XJK * X[ , ll] ) # Same as t(leee) %*% (X[ , jj] * X[ , kk] * X[ , ll])
            coun <- coun + 1
          }
        }
      }
      der3 <- do.call("c", der3)
      
    } # deriv = 3
  }   #         2
  
  out <- list("d1" = der1, "d2" = der2, "d3" = der3)
  
  return( out )
  
}
