##################
#'
#' Log-likelihood of the Tweedie location-scale-shape distribution
#'
#' @description Log-likelihood of the Tweedie location-scale-shape model, and its
#'              derivatives with respect to \code{param = cbind(mu, theta, rho)}, where
#'              \code{mu} is the mean, \code{theta} determines the Tweedie power parameter
#'              via \code{p = (a + b * exp(theta)) / (1 + exp(theta))} and
#'              \code{rho = log(phi)} is the log dispersion. This is a thin wrapper around
#'              \code{mgcv::ldTweedie}, which already returns derivatives w.r.t. \code{mu},
#'              \code{theta} and \code{rho} directly, reshaped into the list format used by
#'              [`gamFactory::llk_gaussian`] and friends. Only 1st and 2nd order
#'              derivatives are available (as in \code{mgcv::twlss}); \code{deriv} greater
#'              than 2 is not supported. See \code{mgcv::ldTweedie} and \code{mgcv::twlss}
#'              for background on the parametrization.
#' @param y a vector of observations (non-negative).
#' @param param a matrix (or list) with 3 columns (elements), containing \code{mu},
#'              \code{theta} and \code{rho}, in this order.
#' @param deriv integer between 0 and 2 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-2
#'              additionally return \code{d1}-\code{d2} (3rd and 4th order derivatives
#'              are not available for this family, see description above).
#' @param a lower limit parameter (> 1) used in the definition of the Tweedie power
#'          parameter \code{p} from \code{theta}. See \code{mgcv::ldTweedie}.
#' @param b upper limit parameter (< 2) used in the definition of the Tweedie power
#'          parameter \code{p} from \code{theta}. See \code{mgcv::ldTweedie}.
#' @rdname llk_twlss
#' @export llk_twlss
#' @examples
#' library(gamFactory)
#'
#' n <- 10
#' param <- cbind(mu = runif(n, 1, 3), theta = rnorm(n, 0, 0.3), rho = rnorm(n, 0, 0.3))
#' y <- rgamma(n, shape = 2, rate = 2)
#'
#' llk_twlss(y = y, param = param, deriv = 2)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_twlss(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_twlss(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_twlss(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 2)$d2))
#'   })
#'
#' check_deriv(obj = obj, param = c(2, 0, 0), ord = 1:2)
#' @export
#'
llk_twlss <- function(y, param, deriv = 0, a = 1.01, b = 1.99, ...) {

  if (is.list(param)) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 3) stop("Wrong number of parameters provided")

  n <- length(y)
  if (nrow(param) == 1 && n > 1) {
    param <- matrix(param, n, 3, byrow = TRUE)
  }

  mu    <- param[ , 1, drop = TRUE]
  theta <- param[ , 2, drop = TRUE]
  rho   <- param[ , 3, drop = TRUE]

  ld <- mgcv::ldTweedie(y, mu = mu, rho = rho, theta = theta, a = a, b = b, all.derivs = TRUE)

  out <- list("d0" = ld[ , 1])

  if (deriv > 0) {
    out[["d1"]] <- list(ld[ , 7], ld[ , 4], ld[ , 2])

    if (deriv > 1) {
      out[["d2"]] <- list(ld[ , 8], ld[ , 9], ld[ , 10], ld[ , 5], ld[ , 6], ld[ , 3])

      if (deriv > 2) {
        stop("llk_twlss only supports derivatives up to order 2 (as in mgcv::twlss)")
      }
    }
  }

  return( out )

}
