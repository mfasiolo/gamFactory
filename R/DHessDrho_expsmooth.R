#'
#' Derivative of exp smooth Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.expsmooth
#' @rdname DHessDrho.expsmooth
#' @export 
#'
DHessDrho.expsmooth <- function(o, llk, DbDr){
  
  d1H <- list()
  
  na <- o$na
  f1 <- o$store$f1; f2 <- o$store$f2; f3 <- o$store$f3
  Xi <- o$store$Xi
  X <- o$store$X0; X1 <- o$store$X1; X2 <- o$store$X2
  le <- llk$d1; lee <- llk$d2; leee <- llk$d3
  
  lggg <- le * f3 + 3 * lee * f1 * f2 + leee * f1^3
  lgg <- le * f2 + lee * f1^2
  leg <- lee * f1
  leeg <- leee * f1
  legg <- leee * f1^2 + lee * f2
  lg <- le * f1
  
  g <- o$store$g; g1 <- o$store$g1; g2 <- o$store$g2; g3 <- o$store$g3
  
  m <- ncol( DbDr )   # Number of smoothing parameters
  
  # Derivatives of alpha and beta wrt log-smoothing parameters
  dArho = DbDr[1:na, , drop = FALSE]
  dBrho = DbDr[-(1:na), , drop = FALSE] 
  
  aaB <- (legg * X + 2 * leg * X1 + le * X2)
  lgggXi <- lggg * g1
  leeeX <- leee * X
  leegXi <- leeg * g1
  leeXi <- lee * g1
  leegX <- leeg * X
  leeX1 <- lee * X1
  leggXi <- legg * g1
  leeX <- lee * X 
  legXi <- leg * g1 
  leeX <- lee * X
  leXi <- le * g1
  lggXi <- lgg * g1
  legX <- leg * X 
  leX1 <- le * X1
  
  .contract3 <- function(X, v, ind){
    d <- length(v)
    Xv <- matrix(nrow = nrow(X), ncol = (d+1)*d/2)
    zz <- 1
    for(jj in 1:d){
      for(kk in jj:d){
      Xv[ , zz] <- X[ , drop(ind[jj, kk, ]), drop = FALSE] %*% v
      zz <- zz + 1
      }
    }
    return( Xv )
  }
  
  .expand2 <- function(X1, X2, ind){
    d <- ncol(v)
    Xv <- matrix(nrow = nrow(X), ncol = d)
    for(ii in 1:d){
      Xv[ , ii] <- X[ , drop(ind[ii, ]), drop = FALSE] %*% v
    }
    return( Xv )
  }
  
  ind  <- trind.generator(na)
  
  for(ii in 1:m){
    dAi <- dArho[ , ii]
    dBi <- dBrho[ , ii]
  
    Q <- .contract2(X = g2, v = dAi, ind = ind$i2)
    Vaaa1_Vaab1 <- drop(lgggXi %*% dAi + aaB %*% dBi)
    UVaaa2Q <- t(g1) %*% (lgg * Q) 
    Vaaa3 <- .vec_to_sym_mat(colSums(drop(lggXi %*% dAi) * g2), na)
    Vaaa4 <- .vec_to_sym_mat(colSums(lg * .contract3(X = g3, v = dAi, ind = ind$i3)), na)
    Vaab2 <- .vec_to_sym_mat(colSums(drop(legX %*% dBi + leX1 %*% dBi) * g2), na)
    derHAA <- t(g1) %*% (Vaaa1_Vaab1 * g1) + UVaaa2Q + t(UVaaa2Q) + Vaaa3 + Vaaa4 + Vaab2
    
    Vbbb1_Vaaa1 <- drop(leeeX %*% dBi + leegXi %*% dAi)
    X1t_Vbba2_X <- t(X1) %*% (drop(leeXi %*% dAi) * X)
    derHBB <- t(X) %*% (Vbbb1_Vaaa1 * X) + X1t_Vbba2_X + t(X1t_Vbba2_X)
    
    Vabb1_Vaba1 <- drop(leegX %*% dBi + leeX1 %*% dBi + leggXi %*% dAi)
    Vabb2_Vaba2 <- drop(leeX %*% dBi + 2 * legXi %*% dAi)
    Vaba3 <- drop(leXi %*% dAi)
    Vaba4 <- leg
    Vaba5 <- le
    derHBA <- t(X) %*% (Vabb1_Vaba1 * g1) + t(X1) %*% (Vabb2_Vaba2 * g1) + 
              t(X2) %*% (Vaba3 * g1) + t(X) %*% (Vaba4 * Q) + t(X1) %*% (Vaba5 * Q)
    
    d1H[[ii]] <- rbind(cbind(derHAA, t(derHBA)), cbind(derHBA, derHBB))
  }
  
  return(d1H)
  
}
