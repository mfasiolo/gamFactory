#'
#' Get B-spline design matrix constructor
#' 
#' @param knots positions of all the knots, both inner and outer.
#' @param m order of the B-spline basis. Here \code{m = 0} for piecewise linear,
#'          \code{m = 1} for quadratic and so on. Number of continuous derivatives
#'          is \code{m}, e.g., if \code{m = 4} fourth order derivatives are continuous
#'          but fifth are not. 
#' 
#' @return A function that take arguments \code{x} and \code{deriv} and returns 
#'         a B-splines model matrix and its first \code{deriv} derivatives w.r.t. \code{x}.
#' 
#' @name basis_bspline
#' @rdname basis_bspline
#' @importFrom splines spline.des
#' @export basis_bspline
#' @examples
#' library(gamFactory)
#' 
#' xseq <- seq(-6, 6, 0.01)
#' k <- 10             
#' m <- 4
#' ko <- k + m + 2
#' knots <- seq(xseq[1], tail(xseq,1), length.out = ko)
#' tmp <- basis_bspline(knots = knots, m = m)$evalX(x = xseq, deriv = 4)
#' 
#' par(mfrow = c(3, 2)) 
#' matplot(xseq, tmp$X0, type = 'l') # Design matrix
#' matplot(xseq, tmp$X1, type = 'l') # 1st derivative w.r.t. x 
#' matplot(xseq, tmp$X2, type = 'l') # 2nd derivative w.r.t. x 
#' matplot(xseq, tmp$X3, type = 'l') # 3rd derivative w.r.t. x
#' matplot(xseq, tmp$X4, type = 'l') # 4rd derivative w.r.t. x
#' rug(knots)
#' 
#' # Checking if derivatives are right by finite differences: it seems so
#' l0 <- basis_bspline(knots = knots, m = m)$evalX(x = xseq, deriv = 4)
#' lm <- basis_bspline(knots = knots, m = m)$evalX(x = xseq-1e-5, deriv = 4)
#' lp <- basis_bspline(knots = knots, m = m)$evalX(x = xseq+1e-5, deriv = 4)
#' par(mfrow = c(2, 2))
#' plot(l0$X1, (lp$X0 - lm$X0) / 2e-5)
#' abline(0, 1, col = 2)
#' plot(l0$X2, (lp$X1 - lm$X1) / 2e-5)
#' abline(0, 1, col = 2)
#' plot(l0$X3, (lp$X2 - lm$X2) / 2e-5)
#' abline(0, 1, col = 2)
#' plot(l0$X4, (lp$X3 - lm$X3) / 2e-5)
#' abline(0, 1, col = 2)
#' 
basis_bspline <- function(knots, m){
  
  force(knots); force(m);
  
  # Numer of total, outer and inner knots
  k <- length(knots)
  ko <- 2 * (m + 1) 
  ki <- k - ko
  
  # Number of interior knots is kint
  if(ki <= 0) stop("Basis dimension too small for b-spline order")
  
  # Range of inner knots
  krange <- range(knots[ko/2 + 1:ki])

  evalX <- function(x, deriv = 0){
    gamFactory:::.basis_bspline(x = x, knots = knots, m = m, deriv = deriv)
  }
  
  out <- structure(list("evalX" = evalX, "krange" = krange), 
                   class = c("B-spline", "basis"))
  
  return(out)
}


##########################
# Internal to construct P-spline design
# 
.basis_bspline <- function(x, knots, m, deriv = 0){
  
  n <- length(x)
  k <- length(knots)
  
  out <- list()
  for(ii in 1:(1+deriv)){
    de <- ii-1
    out[[paste0("X", de)]] <- splines::spline.des(knots, x = x, ord = m + 2, 
                                                  outer.ok = T, derivs = x*0 + de)$design
  }
  
  return( out )
}