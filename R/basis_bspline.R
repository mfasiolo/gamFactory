#'
#' Get B-spline design matrix constructor
#' 
#' @param knots positions of all the knots, both inner and outer.
#' @param m order of the B-spline basis. Here \code{m = 0} for piecewise linear,
#'          \code{m = 1} for quadratic and so on. Number of continuous derivarives
#'          is \code{m}, e.g., if \code{m = 4} fourth order derivatives are continuous
#'          but fifth are not. 
#' 
#' @return A function that take arguments \code{x} and \code{deriv} and returns a P-Spline model matrix and its first
#'         \code{deriv} derivatives w.r.t. \code{x}.
#' 
#' @name basis_bspline
#' @rdname basis_bspline
#' @export basis_bspline
#' @examples
#' library(gamFactory)
#' 
#' xseq <- seq(-6, 6, 0.01)
#' k <- 10
#' m <- 2
#' ko <- k + m + 2
#' knots <- 2 * qt((1:ko)/(ko+1), df = 3)
#' tmp <- basis_bspline(knots = knots, 
#'                      m = m)$evalX(x = xseq, deriv = 3)
#' 
#' par(mfrow = c(2, 2)) 
#' matplot(xseq, tmp$X0, type = 'l') # Design matrix
#' matplot(xseq, tmp$X1, type = 'l') # 1st derivative w.r.t. x 
#' matplot(xseq, tmp$X2, type = 'l') # 2nd derivative w.r.t. x 
#' matplot(xseq, tmp$X3, type = 'l') # 3rd derivative w.r.t. x
#' rug(knots)
#'
basis_bspline <- function(knots, m){
  
  force(knots); force(m); 
  
  # Number of interior knots is kint
  k <- length(knots)
  kint <- k - m - 2
  if(kint <= 0) stop("Basis dimension too small for b-spline order")
  
  out <- list()
  
  evalX <- function(x, deriv = 0){
    gamFactory:::.basis_bspline(x = x, knots = knots, m = m, deriv = deriv)
  }
  
  out <- structure(list("evalX" = evalX), class = c("B-spline", "basis"))
  
  return(out)
}


##########################
# Internal to construct P-spline design
# 
.basis_bspline <- function(x, knots, m, deriv = 0){
  
  n <- length(x)
  k <- length(knots)
  
  # Get full design matrix using also data outside knots: need to call this to get X1, X2 and X3 
  X0 <- splines::spline.des(knots, x = x, ord = m + 2, outer.ok = T)$design

  X1 <- X2 <- X3 <- NULL
  if(deriv > 0){
    X1 <- splines::spline.des(knots, x = x, ord = m + 2, outer.ok = T, derivs = x*0 + 1)$design
    if(deriv > 1){
      X2 <- splines::spline.des(knots, x = x, ord = m + 2, outer.ok = T, derivs = x*0 + 2)$design
    }
    if(deriv > 2){
      X3 <- splines::spline.des(knots, x = x, ord = m + 2, outer.ok = T, derivs = x*0 + 3)$design
    }
  }
  
  return( list("X0" = X0, "X1" = X1, "X2" = X2, "X3" = X3) )
}