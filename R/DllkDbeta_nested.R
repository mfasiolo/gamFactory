#'
#' Derivatives of nested effect
#' 
#' @rdname DllkDbeta.nested
#' @export DllkDbeta.nested
#' @export
#'
DllkDbeta.nested <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv == 0 ){ return( list() ) }
  
  if( is.null(param) ){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  }
  
  # Need to update the object
  if( is.null(o$param) || !identical(param, o$param) || o$deriv < deriv){
    o <- o$eval(param = param, deriv = deriv)
  }
  
  # Extract stuff from object for convenience
  alpha <- param[ 1:o$na ]
  beta <- param[ -(1:o$na) ]
  
  na <- length( alpha )
  nb <- length( beta ) 
  
  # Outer design matrix and its derivatives
  X <- o$store$X0; X1 <- o$store$X1; X2 <- o$store$X2; X3 <- o$store$X3
  
  # Derivatives of outer linear predictor
  f <- o$f; f1 <- o$store$f1; f2 <- o$store$f2; f3 <- o$store$f3
  
  # Inner linear predictor and its derivatives w.r.t. alpha
  g <- o$store$g; g1 <- o$store$g1; g2 <- o$store$g2; g3 <- o$store$g3
  
  der1 <- der2 <- der3 <- der4 <- NULL
  
  le <- llk$d1
  lg <- le * f1
  
  ll_a <- t(g1) %*% lg
  ll_b <- t(X) %*% le
  
  der1 <- c(ll_a, ll_b)
  
  if( deriv >  1){
    lee <- llk$d2
    lgg <- le * f2 + lee * f1^2
    leg <- lee * f1
    
    tmp <- colSums(g2 * lg) 
    ll_aa <- matrix(0, nrow = na, ncol = na)
    ll_aa[lower.tri(ll_aa, diag=TRUE)] <- tmp
    ll_aa <- t(ll_aa)
    ll_aa[lower.tri(ll_aa, diag=TRUE)] <- tmp
    ll_aa <- ll_aa + t(g1) %*% (lgg * g1) 
    ll_bb <- t(X) %*% (lee * X)    # Same as t(X) %*% diag(le2) %*% X
    ll_ba <- t(X) %*% (leg * g1) + t(X1) %*% (le * g1)
    
    der2 <- rbind(cbind(ll_aa, t(ll_ba)), 
                  cbind(ll_ba, ll_bb))
    
    der3 <- list()
    if( deriv >  2  ){
      ina <- trind.generator(na)
      f3 <- drop( X3%*%beta )
      leee <- llk$d3
      lggg <- le * f3 + 3 * lee * f1 * f2 + leee * f1^3
      leeg <- leee * f1
      legg <- leee * f1^2 + lee * f2
      
      coun <- 1
      for(jj in 1:na){      # AXX
        XJ <- g1[ , jj]
        for(kk in jj:na){   # AAX
          XJK <- XJ * g1[ , kk]
          for(ll in kk:na){ # AAA
            der3[[coun]] <- sum(lggg*XJK*g1[,ll]) + 
                            sum(lgg*(g2[,ina$i2[jj, ll]]*g1[,kk]+
                                     g2[,ina$i2[kk, ll]]*g1[,jj] +
                                     g2[,ina$i2[jj, kk]]*g1[,ll])) + 
                            lg %*% g3[,ina$i3[jj, kk, ll]]
            coun <- coun + 1
          }
          for(ll in 1:nb){  # AAB
            der3[[coun]] <- sum((legg*X[,ll]+2*leg*X1[,ll]+le*X2[,ll]) * XJK) + 
                            sum((leg*X[,ll]+le*X1[,ll])*g2[,ina$i2[jj, kk]])
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
