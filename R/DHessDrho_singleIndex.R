#'
#' Derivative of single index Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.singleIndex
#' @rdname DHessDrho.singleIndex
#' @export 
#'
DHessDrho.singleIndex <- function(o, llk, DbDr){
  
  d1H <- list()

  na <- o$na
  f1 <- o$store$f1
  f2 <- o$store$f2
  f3 <- o$store$f3
  Xi <- o$store$Xi
  X <- o$store$X0
  X1 <- o$store$X1
  X2 <- o$store$X2
  le <- llk$d1
  lee <- llk$d2
  leee <- llk$d3
  
  lggg <- le * f3 + 3 * lee * f1 * f2 + leee * f1^3
  lgg <- le * f2 + lee * f1^2
  leg <- lee * f1
  leeg <- leee * f1
  legg <- leee * f1^2 + lee * f2
  
  m <- ncol( DbDr )   # Number of smoothing parameters
  
  # Derivatives of alpha and beta wrt log-smoothing parameters
  dArho = DbDr[1:na, , drop = FALSE]
  dBrho = DbDr[-(1:na), , drop = FALSE]
  
  aaB <- (legg * X + 2 * leg * X1 + le * X2)
  lgggXi <- lggg * Xi
  leeeX <- leee * X
  leegXi <- leeg * Xi
  leeXi <- lee * Xi
  leegX <- leeg * X
  leeX1 <- lee * X1
  leggXi <- legg * Xi
  leeX <- lee * X 
  legXi <- leg * Xi 
  leeX <- lee * X
  leXi <- le * Xi
  
  for(ii in 1:m){
    dAi <- dArho[ , ii]
    dBi <- dBrho[ , ii]
    Vaa <- drop( aaB %*% dBi + lgggXi %*% dAi )
    derHAA <- t(Xi) %*% (Vaa * Xi)
    
    Vbb1 <- drop(leeeX %*% dBi + leegXi %*% dAi)
    Vbb2 <- drop(leeXi %*% dAi)
    derHBB <- t(X) %*% (Vbb1 * X) + t(X1) %*% (Vbb2 * X) + t(X) %*% (Vbb2 * X1)
    
    Vab1 <- drop(leegX %*% dBi + leeX1 %*% dBi + leggXi %*% dAi)
    Vab2 <- drop(leeX %*% dBi + 2 * legXi %*% dAi)
    Vab3 <- drop(leXi %*% dAi)
    derHBA <- t(X) %*% (Vab1 * Xi) + t(X1) %*% (Vab2 * Xi) + t(X2) %*% (Vab3 * Xi)
    
    d1H[[ii]] <- rbind(cbind(derHAA, t(derHBA)), cbind(derHBA, derHBB))
  }

  return(d1H)
  
}
