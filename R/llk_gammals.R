##################
#'
#' Log-likelihood of the Gamma location-scale distribution
#'
#' @description Log-likelihood of the Gamma location-scale model, and its derivatives
#'              with respect to \code{param = cbind(mu, th)}, where \code{mu = log(mean)}
#'              (note: despite the "identity" link name used by \code{mgcv::gammals} for
#'              this parameter, it is internally the log of the mean) and
#'              \code{th = log(phi)} is the log dispersion (\code{phi = 1 / shape}).
#'              Ported directly from the derivative system used internally by
#'              [mgcv::gammals], but re-expressed as derivatives w.r.t. \code{param}
#'              and returned in the list format used by [gamFactory::llk_gaussian]
#'              and friends. See [mgcv::gammals] for background on the
#'              parametrization.
#' @param y a vector of observations (positive).
#' @param param a matrix (or list) with 2 columns (elements), containing \code{mu} and
#'              \code{th}, in this order.
#' @param deriv non-negating positive integer between 0 and 4, indicating the derivative level requires.
#' @rdname llk_gammals
#' @export llk_gammals
#' @examples
#' library(gamFactory)
#'
#' n <- 10
#' param <- cbind(mu = log(runif(n, 1, 3)), th = rnorm(n, 0, 0.3))
#' y <- rgamma(n, shape = exp(-param[ , 2]), scale = exp(param[ , 1] + param[ , 2]))
#'
#' llk_gammals(y = y, param = param, deriv = 4)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_gammals(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_gammals(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_gammals(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 2)$d2))
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_gammals(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 3)$d3))
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_gammals(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = c(0.5, 0), ord = 1:4)
#' @export
#'
llk_gammals <- function(y, param, deriv = 0, ...) {

  if (is.list(param)) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")

  n <- length(y)
  if (nrow(param) == 1 && n > 1) {
    param <- matrix(param, n, 2, byrow = TRUE)
  }

  mu <- param[ , 1, drop = TRUE]
  th <- param[ , 2, drop = TRUE]

  eth <- exp(-th)
  logy <- log(y)
  ethmu <- exp(-th - mu)
  ethmuy <- ethmu * y
  etlymt <- eth * (logy - mu - th)

  d0 <- etlymt - logy - ethmuy - lgamma(eth)

  out <- list("d0" = d0)

  if (deriv > 0) {
    digeth <- digamma(eth)
    d1_mu <- ethmuy - eth
    d1_th <- -etlymt + ethmuy + eth * digeth - eth

    out[["d1"]] <- list(d1_mu, d1_th)

    if (deriv > 1) {
      eth2 <- eth^2
      treth <- trigamma(eth)
      d2_mumu <- -ethmuy
      d2_muth <- eth - ethmuy
      d2_thth <- etlymt - ethmuy - treth * eth2 - eth * digeth + 2 * eth

      out[["d2"]] <- list(d2_mumu, d2_muth, d2_thth)

      if (deriv > 2) {
        eth3 <- eth2 * eth
        g3eth <- psigamma(eth, deriv = 2)
        d3_mumumu <- ethmuy
        d3_mumuth <- ethmuy
        d3_muthth <- ethmuy - eth
        d3_ththth <- -etlymt + ethmuy + g3eth * eth3 + 3 * treth * eth2 + eth * digeth - 3 * eth

        out[["d3"]] <- list(d3_mumumu, d3_mumuth, d3_muthth, d3_ththth)

        if (deriv > 3) {
          eth4 <- eth3 * eth
          d4_mumumumu <- -ethmuy
          d4_mumumuth <- -ethmuy
          d4_mumuthth <- -ethmuy
          d4_muththth <- eth - ethmuy
          d4_thththth <- etlymt - ethmuy - psigamma(eth, deriv = 3) * eth4 - 6 * g3eth * eth3 -
            7 * treth * eth2 - eth * digeth + 4 * eth

          out[["d4"]] <- list(d4_mumumumu, d4_mumumuth, d4_mumuthth, d4_muththth, d4_thththth)
        }
      }
    }
  }

  return( out )

}
