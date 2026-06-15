#' Extract Single-Index Coefficients
#'
#' @description 
#' Extracts the estimated inner single-index coefficients \eqn{\boldsymbol{\alpha}} 
#' from a fitted nested single-index smooth object. The function automatically 
#' accounts for reparametrisation and standardises the direction of the 
#' coefficients by ensuring the first element is always positive.
#'
#' @param object A smooth object containing the single-index information. 
#'   Currently supported classes are \code{"si"} and \code{"si_nexpsm"}. 
#'   This is typically extracted from a fitted model via \code{fit$smooth[[i]]}.
#' @param ... Additional arguments passed to other methods.
#' 
#' @details 
#' Depending on the specific smooth class, the extraction relies on different 
#' internal list elements (e.g., \code{alpha} and \code{a0} for \code{"si"}, 
#' or \code{alpha_si} and \code{alpha_center} for \code{"si_nexpsm"}). 
#' The method seamlessly unifies these differences and returns the final projected coefficients.
#'
#' @return A numeric vector representing the estimated single-index coefficients.
#' 
#' @seealso \code{\link{smooth.construct.si.smooth.spec}}, \code{\link{smooth.construct.si_nexpsm.smooth.spec}}
#'
#' @export
get_alpha <- function(object, ...) {
  UseMethod("get_alpha")
}

#' @rdname get_alpha
#' @export 
get_alpha.si <- function(object, ...) {
  
  si <- object$xt$si
  a0 <- si$a0
  if( is.null(a0) ) a0 <- si$alpha * 0
  
  a <- si$B %*% (si$alpha + a0)
  drop(sign(a[1]) * a)
}

#' @rdname get_alpha
#' @export 
get_alpha.si_nexpsm <- function(object, ...) {
  
  si <- object$xt$si
  alpha_center <- si$alpha_center
  if( is.null(alpha_center) ) alpha_center <- si$alpha_si * 0
  
  a <- si$B_si %*% (si$alpha_si + alpha_center)
  drop(sign(a[1]) * a)
}