.vec_to_sym_mat <- function(x, d){
  if(!is.vector(x)){ x <- drop(x) }
  m1 <- matrix(NA, d, d)
  m1[ lower.tri(m1, diag=TRUE) ] <- x
  m1 <- pmax(m1, t(m1), na.rm=TRUE)
  return(m1)
}