#' Build Nested Adaptive Exponential Smooth Effect Evaluator
#' 
#' @description This function acts as a constructor for nested adaptive 
#' exponential smoothing. It bundles the target sequence, smoothing covariates, initial state, 
#' and an outer basis expansion method into an object that 
#' can evaluate the effect and its analytical derivatives dynamically.
#' 
#' @param y A numeric vector representing the time series or sequence data to be 
#' exponentially smoothed.
#' @param Xi A numeric design matrix of dimensions \code{n} by \code{p} containing 
#' the historical covariates that drive the time-varying smoothing rate.
#' @param basis A structure or list containing an outer basis evaluation function 
#' \code{evalX(x, deriv)} which evaluates the final spline transformation stage.
#' @param x0 A numeric scalar representing the initial state of the smooth at 
#' time step 0. Defaults to \code{NULL}, which sets the initial state to \code{y[1]}.
#' @param times An optional integer vector of indices specifying the precise time 
#' points at which to return evaluations and analytical derivatives.
#' 
#' @details The function maps a two-layer hierarchical smoothing architecture:
#' \enumerate{
#'   \item \strong{Inner Adaptive Exponential Smooth:} The raw parameters are split. The first 
#'         \code{ncol(Xi) + 1} parameters (\code{alpha}) are fed directly into 
#'         \code{\link{expsmooth}} to yield a time-varying smoothed trajectory vector:
#'         \deqn{\eta_i = \text{expsmooth}(y, \mathbf{X}_i, \bm{\alpha})_i}
#'   \item \strong{Outer Spline Expansion:} The smooth vector elements \eqn{\eta} 
#'         are evaluated through the outer \code{basis} transformation at the observation points
#'         corresponding to \code{times}. The 
#'         remaining parameter elements (\code{beta}) serve as linear weights 
#'         against this expanded basis matrix:
#'         \deqn{f(\eta_i) = \mathbf{B}_0(\eta_i) \bm{\beta}}
#' }
#' @return An object of class \code{c("nexpsm", "nested")} containing a locked-environment closure:
#' \itemize{
#'   \item{\code{.eval(param, deriv = 0)}}{ An evaluation function that takes a combined 
#'         numeric parameter vector \code{param} (where the first \code{ncol(Xi) + 1} 
#'         elements represent the smoothing parameters, and the remainder represent 
#'         the outer spline coefficients) and returns an updated effect state evaluation.}
#' }
#' @seealso [trans_nexpsm], [s_nest]
#' @name eff_nexpsm
#' @rdname eff_nexpsm
#' @export eff_nexpsm
#'
eff_nexpsm <- function(y, Xi, basis, x0 = NULL, times = NULL){
  
  force(y); force(Xi); force(basis); force(x0); force(times);
  
  na <- ncol(Xi) + 1
  
  incall <- as.expression(quote(do.call("expsmooth", list("y" = y, "Xi" = Xi, "beta" = alpha, "times" = times, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_nexpsm", list("y" = y, "Xi" = Xi, "basis" = basis, "x0" = x0, "times" = times), quote = TRUE)))
  
  .eval <- .get_eff_eval_general()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("nexpsm", "nested"))
  
  return( out )
  
}









