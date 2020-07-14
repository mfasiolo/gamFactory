
.Hess.standard_standard <- function(o1, o2, llk){
  
  out <- crossprod(o2$store$X, llk$d2 * o1$store$X)
  
  return( out )
  
}


.Hess.singleIndex_standard <- function(o1, o2, llk){

  out <- cbind(crossprod(o2$store$X, llk$d2 * o1$store$f1 * o1$store$Xi), 
               crossprod(o2$store$X, llk$d2 * o1$store$X0))
    
 return( out )
  
}

.Hess.standard_singleIndex <- function(o1, o2, llk){

  # Swapped o1 with o2 and lp1 with lp2
  t( .Hess.singleIndex_standard(o2, o1, llk) ) 
  
}

.Hess.singleIndex_singleIndex <- function(o1, o2, llk){
  
  lgg1 <- llk$d2 * o1$store$f1 * o2$store$f1
  
  # o2 vertical, o1 horizontal
  out <- rbind(cbind(crossprod(o2$store$Xi, lgg1 * o1$store$Xi), 
                     crossprod(o2$store$Xi, llk$d2 * o2$store$f1 * o1$store$X0)),
               cbind(crossprod(o2$store$X0, (llk$d2 * o1$store$f1) * o1$store$Xi), 
                     crossprod(o2$store$X0, llk$d2 * o1$store$X0)))
  
  
  return( out )
  
}