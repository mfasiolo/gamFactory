#'
#' Specifying a MGKS transformation
#' 
#' @name trans_mgks
#' @rdname trans_xxx
#' @export trans_mgks
#'
trans_mgks <- function(n0, alpha){
 
  if( missing(n0) ){ stop("Argument \"n0\" is missing but a value is required") }
   
  out <- as.list(match.call())[-1]
  out$type <- "mgks"
 
  return(out)
  
} 

#' @rdname trans_xxx
#' @export trans_nexpsm
#'
trans_nexpsm <- function(S, alpha){
  
  out <- as.list(match.call())[-1]
  out$type <- "nexpsm"
  
  return(out)
  
} 

#' @rdname trans_xxx
#' @export trans_linear
#'
trans_linear <- function(pord, S, alpha){
  
  out <- as.list(match.call())[-1]
  out$type <- "si"
  
  return(out)
  
}