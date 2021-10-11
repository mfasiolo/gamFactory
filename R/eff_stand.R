#'
#' Build standard effects
#' 
#' @param X model matrix, such that effect will be \code{f = X \%*\% beta}.
#' @name eff_stand
#' @rdname eff_stand
#' @export eff_stand
#'
eff_stand <- function(X){
  
  force(X)
  
  eval <- function(param, deriv = 0){
    
    o <- eff_stand(X = X)
    o$f <- drop( X %*% param ) 
    if(deriv >= 1){
      o$store <- list("X" = X)
    }
    o$deriv <- deriv
    o$param <- param
  
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("stand"))
  
}

