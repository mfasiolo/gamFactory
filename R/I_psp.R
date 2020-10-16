# P-spline penalty of order "ord" and dimension "d"
.psp <- function(d, ord){
  S <- diag(1, d)
  if( ord > 0 ){
    S <- diff(S, differences = ord)
    S <- crossprod(S)
  }
  return(S)
}