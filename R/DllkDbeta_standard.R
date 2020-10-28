#'
#' Derivatives of log-lik w.r.t. coeffs of a standard effect
#' 
#' @rdname DllkDbeta.standard
#' @export DllkDbeta.standard
#' @export
#'
DllkDbeta.standard <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv == 0 ){ return( list() ) }
  
  X <- o$store$X
  n <- nrow(X)
  p <- ncol( X )
  
  d1 <- colSums(llk$d1 * X)
  
  d2 <- d3 <- NULL
  if( deriv >= 2){
    d2 <- crossprod(X, llk$d2 * X)
  if( deriv >=  3 ){
    d3 <- list()
    coun <- 1
    for(jj in 1:p){
      XJ <- llk$d3 * X[ , jj]
      for(kk in jj:p){
        XJK <- XJ * X[ , kk]
        for(ll in kk:p){
          # Same as t(leee) %*% (X[ , jj] * X[ , kk] * X[ , ll])
          d3[[coun]] <- sum( XJK * X[ , ll] ) 
          coun <- coun + 1
        }
      }
    }
    d3 <- do.call("c", d3)
  } # deriv >= 3
  } # deriv >= 2
  
  
  out <- list("d1" = d1, "d2" = d2, "d3" = d3)
  
  return( out )
  
}