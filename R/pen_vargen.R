s#'
#' Penalty on variance of a non-linear effect
#' 
#' @description Penalty on variance of a single index vector. 
#' @param a single index parameter vector of dimension p.
#' @param x an n by p matrix of covariates to be projected.
#' @param v target variance of the projected covariates.
#' @param deriv number of derivatives to be calculated.
#' @param DaDr derivative for single index vector w.r.t. smoothing parameter.
#' @name pen_vargen
#' @rdname pen_vargen
#' @export pen_vargen
#'
pen_vargen <- function(o, v, deriv = 0){
  
  p <- o$na
  
  stopifnot( !is.null(o$param[1:p]) || o$deriv < deriv )
  
  # Need to update the object
  g <- o$store$g
  n <- length(g)
  
  # Empirical variance
  vhat <- sum(g^2)/n - mean(g)^2
  
  # Loss
  l0 <- (vhat - v)^2 
  
  l1 <- l2 <- l3 <- NULL
  if( deriv ){
    g1 <- o$store$g1
    g1cT <- t(g1) - colMeans(g1)
    gc <- g - mean(g)
    q <- g1cT %*% gc
    l1 <- 4/n * (vhat - v) * q
    if( deriv > 1 ){
      g2cT <- t(o$store$g2) - colMeans(o$store$g2)
      CP <- tcrossprod(g1cT) + .vec_to_sym_mat(g2cT %*% gc, p)
      l2 <- 4/n * (2/n * tcrossprod(q) + (vhat - v) * CP) 
      if( deriv > 2 ){
        ind  <- trind.generator(p)
        A <- B <- numeric(p^2 * (p + 1) / 2)
        l3 <- numeric(choose(p + 2, p-1))
        g3cT <- t(o$store$g3) - colMeans(o$store$g3)
        g3cTgc <- g3cT %*% gc
        kk <- 1
        for(i1 in 1:p){
          for(i2 in 1:p){ # Yes 1:p, NOT i1:p !!
            for(i3 in i2:p){
              A[kk] <- q[i1]*CP[i2, i3]
              B[kk] <- crossprod(g1cT[i1, ], g2cT[ind$i2[i2, i3], ])
              kk <- kk + 1
            }
          }
        }
        off <- p * (p + 1) / 2
        mid <- function(a,b,c){ return( off*(a-1) + ind$i2[b, c] ) } 
        kk <- 1
        for(i1 in 1:p){
          for(i2 in i1:p){ 
            for(i3 in i2:p){
              l3[kk] <- 8/n^2 * (A[mid(i1,i2,i3)]+A[mid(i2,i1,i3)]+A[mid(i3,i1,i2)]) + 
                        4/n * (vhat - v) * (B[mid(i1,i2,i3)]+B[mid(i2,i1,i3)]+B[mid(i3,i1,i2)]+g3cTgc[kk])
              kk <- kk + 1
            }
          }
        }
      }
    }
  } 
  
  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = l3) )
  
}


#' @name pen_vargen
#' @rdname pen_vargen
#' @export pen_vargen_outer
pen_vargen_outer <- function(o, v, DaDr){
  
  g <- o$store$g
  g1 <- o$store$g1
  g2 <- o$store$g2
  g3 <- o$store$g3
  m <- ncol(DaDr)
  n <- nrow(g1)  
  p <- ncol(g1)
  
  vhat <- sum(g^2)/n - mean(g)^2
  g1cT <- t(g1) - colMeans(g1)
  gc <- g - mean(g)
  q <- g1cT %*% gc
  g2cT <- t(g2) - colMeans(g2)
  g2c <- t(g2cT)
  g3cvett <- (t(g3) - colMeans(g3)) %*% gc
  CP <- tcrossprod(g1cT) + .vec_to_sym_mat(g2cT %*% gc, p)
  
  ind  <- trind.generator(p)
  d1H <- list()
  for(ii in 1:m) {
    dar <- DaDr[ , ii]
    CPa <- CP %*% dar
    Z <- .contract2(X = g2c, v = dar, ind = ind$i2) 
    Zg1c <- g1cT %*% Z
    d1H[[ii]] <- 8 / n^2 * ( tcrossprod(CPa, q) + tcrossprod(q, CPa) + CP* drop(crossprod(dar, q)) ) + 
                 4 / n * (vhat - v) *  ( Zg1c + t(Zg1c) + 
                                         .vec_to_sym_mat(g2cT %*% crossprod(g1cT, dar), p) + 
                                         .vec_to_sym_mat(.contract3_vect(g3cvett, dar, ind$i3), p) )
  } 
  
  return( d1H ) 
  
}