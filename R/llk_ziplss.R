##################
#'
#' Log-likelihood of the Zero-Inflated Poisson (ZIP) distribution
#'
#' @description Log-likelihood of the zero-inflated Poisson location-scale model, and its
#'              derivatives with respect to \code{param = cbind(g, eta)}, where \code{g}
#'              is the linear predictor for \code{log(lambda)} (\code{lambda} being the
#'              Poisson rate of the "at risk" component) and \code{eta} is the linear
#'              predictor for the zero-inflation probability, related via
#'              \code{p = 1 - exp(-exp(eta))} (a complementary log-log style link). This
#'              is simply a thin wrapper around \code{mgcv::zipll}, which already returns
#'              derivatives w.r.t. \code{g} and \code{eta} directly, reshaped into the list
#'              format used by [`gamFactory::llk_gaussian`] and friends. See
#'              \code{mgcv::ziplss} for background on the parametrization.
#' @param y a vector of observations (non-negative integers).
#' @param param a matrix (or list) with 2 columns (elements), containing \code{g} and
#'              \code{eta}, in this order.
#' @param deriv integer between 0 and 4 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-4
#'              additionally return \code{d1}-\code{d4}.
#' @rdname llk_ziplss
#' @export llk_ziplss
#' @examples
#' library(gamFactory)
#'
#' n <- 10
#' param <- cbind(g = rnorm(n, 1, 0.3), eta = rnorm(n, 0.5, 0.3))
#' lambda <- exp(param[ , 1])
#' p <- 1 - exp(-exp(param[ , 2]))
#' y <- ifelse(p > runif(n), rpois(n, lambda), 0)
#'
#' llk_ziplss(y = y, param = param, deriv = 4)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_ziplss(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_ziplss(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_ziplss(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 2)$d2))
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_ziplss(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 3)$d3))
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_ziplss(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = c(1, 0.5), ord = 1:4)
#' @export
#'
llk_ziplss <- function(y, param, deriv = 0, ...) {

  if (is.list(param)) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")

  n <- length(y)
  if (nrow(param) == 1 && n > 1) {
    param <- matrix(param, n, 2, byrow = TRUE)
  }

  g   <- param[ , 1, drop = TRUE]
  eta <- param[ , 2, drop = TRUE]

  zl <- mgcv::zipll(y, g, eta, deriv = deriv)

  out <- list("d0" = zl$l)

  if (deriv > 0) {
    out[["d1"]] <- list(zl$l1[ , 1], zl$l1[ , 2])

    if (deriv > 1) {
      out[["d2"]] <- list(zl$l2[ , 1], zl$l2[ , 2], zl$l2[ , 3])

      if (deriv > 2) {
        out[["d3"]] <- list(zl$l3[ , 1], zl$l3[ , 2], zl$l3[ , 3], zl$l3[ , 4])

        if (deriv > 3) {
          out[["d4"]] <- list(zl$l4[ , 1], zl$l4[ , 2], zl$l4[ , 3], zl$l4[ , 4], zl$l4[ , 5])
        }
      }
    }
  }

  return( out )

}
