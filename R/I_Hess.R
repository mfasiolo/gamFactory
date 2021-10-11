
.Hess.stand_stand <- function(o1, o2, llk){
  
  out <- crossprod(o2$store$X, llk$d2 * o1$store$X)
  
  return( out )
  
}


.Hess.si_stand <- .Hess.nexpsm_stand <- function(o1, o2, llk){

  out <- cbind(crossprod(o2$store$X, llk$d2 * o1$store$f1 * o1$store$g1), 
               crossprod(o2$store$X, llk$d2 * o1$store$X0))
    
 return( out )
  
}

.Hess.stand_si <- .Hess.stand_nexpsm <- function(o1, o2, llk){

  # Swapped o1 with o2 and lp1 with lp2
  t( .Hess.si_stand(o2, o1, llk) ) 
  
}

.Hess.si_nexpsm <- function(o1, o2, llk){
  
  lgg1 <- llk$d2 * o1$store$f1 * o2$store$f1
  
  # o2 vertical, o1 horizontal
  out <- rbind(cbind(crossprod(o2$store$g1, lgg1 * o1$store$g1), 
                     crossprod(o2$store$g1, llk$d2 * o2$store$f1 * o1$store$X0)),
               cbind(crossprod(o2$store$X0, (llk$d2 * o1$store$f1) * o1$store$g1), 
                     crossprod(o2$store$X0, llk$d2 * o1$store$X0)))
  
  
  return( out )
  
}

.Hess.nexpsm_si <- function(o1, o2, llk){
  
  # Swapped o1 with o2 and lp1 with lp2
  out <- t( .Hess.si_nexpsm(o1 = o2, o2 = o1, llk) )
  
  return( out )
  
}

.Hess.si_si <- .Hess.nexpsm_nexpsm <- function(o1, o2, llk){
  
  lgg1 <- llk$d2 * o1$store$f1 * o2$store$f1
  
  # o2 vertical, o1 horizontal
  out <- rbind(cbind(crossprod(o2$store$g1, lgg1 * o1$store$g1), 
                     crossprod(o2$store$g1, llk$d2 * o2$store$f1 * o1$store$X0)),
               cbind(crossprod(o2$store$X0, (llk$d2 * o1$store$f1) * o1$store$g1), 
                     crossprod(o2$store$X0, llk$d2 * o1$store$X0)))
  
  
  return( out )
  
}

