#'
#' Transform symmetric matrix to vector 
#' 
#' @description The input is a symmetric matrix, the output is a vector of its unique elements by row. 
#' @name symmat2vec
#' @param H a symmetric matrix.
#' @return A vector of unique elements of H, by row.
#' @rdname symmat2vec
#' @export symmat2vec
#' @examples
#' 1:3 %*% t(1:3)
#' symmat2vec(1:3 %*% t(1:3))
#' 
symmat2vec <- function(H){
  d <- nrow( H )
  h <- numeric(d*(d+1)/2)
  kk <- 1
  for(ir in 1:d){
    for(ic in ir:d){
    h[kk] <- H[ir, ic] 
    kk <- kk + 1
    }
  } 
  return(  h )
}