
.Hess.singleIndex_standard <- function(o1, o2, llk){

  out <- cbind(crossprod(o2$xtra$X, llk$d2 * o1$xtra$f1 * o1$xtra$Xi), 
               crossprod(o2$xtra$X, llk$d2 * o1$xtra$X))
    
 return( out )
  
}

.Hess.standard_singleIndex <- function(o1, o2, llk){

  # Swapped o1 with o2 and lp1 with lp2
  t( .Hess.singleIndex_standard(o2, o1, llk) ) 
  
}

.Hess.singleIndex_singleIndex <- function(o1, o2, llk){
  
  lgg1 <- llk$d2 * o1$xtra$f1 * o2$xtra$f1
  
  # o2 vertical, o1 horizontal
  out <- rbind(cbind(crossprod(o2$xtra$Xi, lgg1 * o1$xtra$Xi), 
                     crossprod(o2$xtra$Xi, (llk$d2 * o2$xtra$f1) * o1$xtra$X)),
               cbind(crossprod(o2$xtra$X, (llk$d2 * o1$xtra$f1) * o1$xtra$Xi), 
                     crossprod(o2$xtra$X, llk$d2 * o1$xtra$X)))
  
  
  return( out )
  
}