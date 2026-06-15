#' Build Standard Effects
#' 
#' @description This function acts as a constructor for standard parametric linear 
#' or fixed-basis smooth effects. 
#' 
#' @param X A numeric model matrix of dimensions \code{n} by \code{p}, where rows represent 
#' observations and columns represent parametric covariates or fixed spline basis functions.
#' 
#' @details The function maps a conventional linear model component.
#' For a given regression coefficient vector \code{param} (\eqn{\bm{\beta}}):
#' \deqn{f = \mathbf{X} \bm{\beta}}
#' 
#' When \code{deriv >= 1}, the matrix \code{X} is packed directly into the \code{store} 
#' slot. Because the effect is purely linear with respect to \code{param}, its first-order 
#' analytical derivative matrix is identically \code{X}, and all higher-order derivatives 
#' are zero.
#' 
#' @return An object of class \code{"stand"} containing a single environment-locked closure:
#' \itemize{
#'   \item{\code{eval(param, deriv = 0)}}{ An evaluation function that takes a numeric 
#'         coefficient vector \code{param} and returns a calculated state object containing 
#'         the evaluated linear predictor payload and design matrix store.}
#' }
#' 
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
  
  return(out)
  
}

