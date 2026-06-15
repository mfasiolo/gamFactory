#' Build Multivariate Kernel Smooth Effect Evaluator
#' 
#' @description This function acts as a constructor for nested multivariate 
#' kernel smooth (MGKS) effects. It bundles the target response sequence, 
#' pre-computed spatial distance matrices, and an outer basis expansion method
#' into an object that can evaluate the combined effect and its higher-order analytical derivatives.
#' 
#' @param y A numeric vector or matrix representing the response variable baseline 
#' data passed to the kernel smoother.
#' @param dist A list of coordinate-wise squared distance matrices across the 
#' target covariate dimensions used to compute the Gaussian kernel weights.
#' @param basis A structure or list containing an outer basis evaluation function 
#' \code{evalX(x, deriv)} which evaluates the final spline transformation stage.
#' 
#' @details The function establishes a two-layer nested effect structure:
#' \enumerate{
#'   \item \strong{Inner Transformation:} The parameter vector is partitioned. 
#'         The first \code{na} parameters (\code{alpha}) define the log-inverse 
#'         bandwidths used by the underlying \code{mgks} function to compute an 
#'         adaptive multivariate kernel smooth baseline:
#'         \deqn{\tilde{s}(\mathbf{x}_i) = \text{mgks}(y, \text{dist},\alpha)_i.}
#'   \item \strong{Outer Smooth Expansion:} The calculated smooth values \eqn{\eta} 
#'         are dynamically fed into the outer \code{basis} object. The remaining 
#'         parameters (\code{beta}) are used as linear coefficients against this 
#'         basis expansion matrix:
#'         \deqn{f(\tilde{s}(\mathbf{x}_i) ) = \mathbf{B}_0(\tilde{s}(\mathbf{x}_i) ) \bm{\beta}}
#' }
#' 
#' This composition allows the package's optimization engines to jointly fit 
#' both the non-linear bandwidth profiles and the outer spline curves simultaneously.
#' 
#' @return An object of class \code{c("mgks", "nested")} containing a locked-environment closure:
#' \itemize{
#'   \item{\code{.eval(param, deriv = 0)}}{ An evaluation function that takes a combined 
#'         numeric parameter vector \code{param} (where the first \code{length(dist) + 1} 
#'         elements represent the kernel tuning metrics, and the remainder represent 
#'         the outer spline coefficients) and returns an updated effect state evaluation.}
#' }
#' @seealso [trans_mgks], [s_nest]
#' @name eff_mgks
#' @rdname eff_mgks
#' @export eff_mgks
#'
eff_mgks <- function(y, dist, basis){
  
  force(y); force(dist); force(basis);
  
  na <- length(dist) + 1
  Xi <- NULL
  
  incall <- as.expression(quote(do.call("mgks", list("y" = y, "dist" = dist, "beta" = alpha, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_mgks", list("y" = y, "dist" = dist, "basis" = basis), quote = TRUE)))
  
  .eval <- .get_eff_eval_general()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("mgks", "nested"))
  
  return( out )
  
}

