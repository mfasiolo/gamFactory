##################
#'
#' Log-likelihood of the Gumbel location-scale distribution
#'
#' @description Log-likelihood of the Gumbel location-scale model, and its derivatives
#'              with respect to \code{param = cbind(mu, bt)}, where \code{mu} is the mode
#'              and \code{bt = log(beta)} is the log scale of the Gumbel distribution.
#'              Ported directly from the derivative system used internally by
#'              \code{mgcv::gumbls}, but re-expressed as derivatives w.r.t. \code{param}
#'              and returned in the list format used by [`gamFactory::llk_gaussian`]
#'              and friends. See \code{mgcv::gumbls} for background on the
#'              parametrization.
#' @param y a vector of observations.
#' @param param a matrix (or list) with 2 columns (elements), containing \code{mu} and
#'              \code{bt}, in this order.
#' @param deriv integer between 0 and 4 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-4
#'              additionally return \code{d1}-\code{d4}.
#' @rdname llk_gumbls
#' @export llk_gumbls
#' @examples
#' library(gamFactory)
#'
#' n <- 10
#' param <- cbind(mu = rnorm(n, 2, 0.5), bt = rnorm(n, 0, 0.3))
#' u <- runif(n)
#' y <- param[ , 1] - exp(param[ , 2]) * log(-log(u))
#'
#' llk_gumbls(y = y, param = param, deriv = 4)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_gumbls(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_gumbls(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_gumbls(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 2)$d2))
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_gumbls(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 3)$d3))
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_gumbls(y = y, param = matrix(param, n, 2, byrow = TRUE), deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = c(2, 0), ord = 1:4)
#' @export
#'
llk_gumbls <- function(y, param, deriv = 0, ...) {

  if (is.list(param)) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")

  n <- length(y)
  if (nrow(param) == 1 && n > 1) {
    param <- matrix(param, n, 2, byrow = TRUE)
  }

  mu <- param[ , 1, drop = TRUE]
  bt <- param[ , 2, drop = TRUE]

  eb <- exp(-bt)
  z <- (y - mu) * eb
  ez <- exp(-z)

  d0 <- -bt - z - ez

  out <- list("d0" = d0)

  if (deriv > 0) {
    lz <- ez - 1
    zm <- -eb
    zb <- -z

    d1_mu <- lz * zm
    d1_bt <- lz * zb - 1

    out[["d1"]] <- list(d1_mu, d1_bt)

    if (deriv > 1) {
      lzz <- -ez
      zmb <- eb
      zbb <- z

      d2_mumu <- lzz * zm^2
      d2_mubt <- lzz * zm * zb + lz * zmb
      d2_btbt <- lzz * zb^2 + lz * zbb

      out[["d2"]] <- list(d2_mumu, d2_mubt, d2_btbt)

      if (deriv > 2) {
        lzzz <- ez
        zbbb <- -z
        zmbb <- -eb

        d3_mumumu <- lzzz * zm^3
        d3_mumubt <- lzzz * zm^2 * zb + 2 * lzz * zm * zmb
        d3_mubtbt <- lzzz * zb^2 * zm + 2 * lzz * zb * zmb + lzz * zbb * zm + lz * zmbb
        d3_btbtbt <- lzzz * zb^3 + 3 * lzz * zb * zbb + lz * zbbb

        out[["d3"]] <- list(d3_mumumu, d3_mumubt, d3_mubtbt, d3_btbtbt)

        if (deriv > 3) {
          lzzzz <- -ez
          zbbbb <- z
          zmbbb <- eb

          d4_mumumumu <- lzzzz * zm^4
          d4_mumumubt <- lzzzz * zm^3 * zb + 3 * lzzz * zm^2 * zmb
          d4_mumubtbt <- lzzzz * zm^2 * zb^2 + 4 * lzzz * zm * zb * zmb + lzzz * zm^2 * zbb +
            2 * lzz * zmb^2 + 2 * lzz * zm * zmbb
          d4_mubtbtbt <- lzzzz * zb^3 * zm + 3 * lzzz * zb^2 * zmb + 3 * lzzz * zm * zb * zbb +
            3 * lzz * zmb * zbb + 3 * lzz * zb * zmbb + lzz * zm * zbbb + lz * zmbbb
          d4_btbtbtbt <- lzzzz * zb^4 + 6 * lzzz * zb^2 * zbb + 3 * lzz * zbb^2 + 4 * lzz * zb * zbbb + lz * zbbbb

          out[["d4"]] <- list(d4_mumumumu, d4_mumumubt, d4_mumubtbt, d4_mubtbtbt, d4_btbtbtbt)
        }
      }
    }
  }

  return( out )

}
