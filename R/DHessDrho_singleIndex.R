#'
#' Derivative of single index Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.singleIndex
#' @rdname DHessDrho.singleIndex
#' @export DHessDrho.singleIndex
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
  
  for(ii in 1:m){
    dAi <- dArho[ , ii]
    dBi <- dBrho[ , ii]
    Vaa <- drop( (legg * X + 2 * leg * X1 + le * X2) %*% dBi + lggg * Xi %*% dAi )
    derHAA <- t(Xi) %*% (Vaa * Xi)
    
    Vbb1 <- drop(leee * X %*% dBi + leeg * Xi %*% dAi)
    Vbb2 <- drop(lee * Xi %*% dAi)
    derHBB <- t(X) %*% (Vbb1 * X) + t(X1) %*% (Vbb2 * X) + t(X) %*% (Vbb2 * X1)
    
    Vab1 <- drop(leeg * X %*% dBi + lee * X1 %*% dBi + legg * Xi %*% dAi)
    Vab2 <- drop(lee * X %*% dBi + 2 * leg * Xi %*% dAi)
    Vab3 <- drop(le * Xi %*% dAi)
    derHBA <- t(X) %*% (Vab1 * Xi) + t(X1) %*% (Vab2 * Xi) + t(X2) %*% (Vab3 * Xi)
    
    d1H[[ii]] <- rbind(cbind(derHAA, t(derHBA)), cbind(derHBA, derHBB))
  }

  return(d1H)
  
}
