#'
#' Build single index effect
#' 
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param basis function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @name eff_si
#' @rdname eff_si
#' @export eff_si
#'
eff_si <- function(Xi, basis, a0 = NULL, positive_si = FALSE){
  positive_si <- isTRUE(positive_si)
  force(Xi); force(basis); force(a0); force(positive_si)

  eval <- function(param, deriv = 0){
    
    na <- ncol( Xi )
    nb <- length(param) - na 
    
    # Single index and spline coefficients
    alpha <- param[ 1:na ]
    beta <- param[ -(1:na) ]
    
    if( is.null(a0) ){
      a0 <- alpha * 0
    }

    # Project covariates on single index vector 
    if (positive_si) {
      ax <- drop( Xi %*% exp(alpha + a0) )
    } else {
      ax <- drop( Xi %*% (alpha + a0) )
    }
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- basis$evalX(x = ax, deriv = deriv)
    store$Xi <- Xi
    
    if (positive_si) {
      store$g <- ax
    } 
    
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      
      if (positive_si) {
        g_diag <- t(t(Xi) * as.vector(exp(alpha + a0)))
        store$g1 <- g_diag
      } else {
        store$g1 <- Xi
      }
      
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        
        if (positive_si) {
          # expand to n x [na*(na+1)/2], use 0 to fill the position
          n_lower <- na * (na + 1) / 2
          g2_full <- matrix(0, nrow = nrow(Xi), ncol = n_lower)
          
          idx_mat <- matrix(0, na, na)
          idx_mat[lower.tri(idx_mat, diag = TRUE)] <- 1:n_lower
          diag_indices <- diag(idx_mat)
          
          g2_full[, diag_indices] <- g_diag
          store$g2 <- g2_full
        } 
        
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
          
          if (positive_si) {
            # expand to match required size from general case
            n_tensor <- na * (na + 1) * (na + 2) / 6
            g3_full <- matrix(0, nrow = nrow(Xi), ncol = n_tensor)
            
            diag_indices_g3 <- numeric(na)
            count <- 1
            for(i in 1:na) {
              for(j in i:na) {
                for(k in j:na) {
                  if(i == j && j == k) {
                    diag_indices_g3[i] <- count
                  }
                  count <- count + 1
                }
              }
            }
            
            g3_full[, diag_indices_g3] <- g_diag
            store$g3 <- g3_full
          }
        }
      }
    }
    
    o <- eff_si(Xi = Xi, basis = basis, a0 = a0, positive_si = positive_si)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$a0 <- a0
    o$na <- na
    o$store <- store
    o$deriv <- deriv

    return( o )
    
  }
  
if (positive_si) {
  out <- structure(list("eval" = eval), class = c("si_posi", "nested"))
} else {
  out <- structure(list("eval" = eval), class = c("si", "nested"))
}
  return( out )
  
}









