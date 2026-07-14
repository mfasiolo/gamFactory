.my_eigenfix <- function(H, tol){
  
  eg <- eigen(H)
  
  bad <- eg$values < eg$values[1] * tol
  if( any(bad) ){
    eg$values[bad] <- eg$values[1] * tol
    H <- eg$vectors %*% (eg$values * t(eg$vectors))
  }
  
  return(H)
}