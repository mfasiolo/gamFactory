#'
#' Build standard effects
#' 
#' @name buildStandardEffect
#' @rdname buildStandardEffect
#' @export buildStandardEffect
#'
buildStandardEffect <- function(X){
  
  eval <- function(param, deriv){
    
    o <- buildStandardEffect(X = X)
    o$f <- drop( X %*% param ) 
    if(deriv >= 1){
      o$store <- list("X" = X)
    }
    o$deriv <- deriv
  
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("standard"))
  
}

