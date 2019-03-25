#'
#' Transform symmetric matrix to vector 
#' 
#' @description The input is a symmetric matrix, the output is a vector of its unique elements by row. 
#' @name flattenSymMat
#' @param H a symmetric matrix.
#' @return A vector of unique elements of H, by row.
#' @rdname flattenSymMat
#' @export flattenSymMat
#' 
flattenSymMat <- function(H){
  d <- nrow( H )
  h <- list( )
  for(ir in 1:d){
    h[[ir]] <- drop( H[ir, ir:d] )
  } 
  return( do.call("c", h) )
}