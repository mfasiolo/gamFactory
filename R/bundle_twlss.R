#'
#' Bundle for the Tweedie location-scale-shape model
#'
#' @description Ported from [`mgcv::twlss`]. The three linear predictors correspond,
#'              in order, to the mean \code{mu}, \code{theta} (which determines the
#'              Tweedie power parameter via \code{p = (a + b * exp(theta)) / (1 + exp(theta))})
#'              and \code{rho = log(phi)}, the log dispersion. Only 1st and 2nd order
#'              derivatives are available (as in \code{mgcv::twlss}, whose
#'              \code{available.derivs} is 0).
#' @param a lower limit parameter (> 1) used in the definition of the Tweedie power
#'          parameter \code{p} from \code{theta}. See \code{mgcv::ldTweedie}.
#' @param b upper limit parameter (< 2) used in the definition of the Tweedie power
#'          parameter \code{p} from \code{theta}. See \code{mgcv::ldTweedie}.
#' @name bundle_twlss
#' @rdname bundle_twlss
#' @export
#'
bundle_twlss <- function(a = 1.01, b = 1.99){
  force(a)
  force(b)

  out <- list(np = 3,
              available_deriv = 2,
              llk = function(y, param, deriv = 0, ...) {
                llk_twlss(y = y, param = param, deriv = deriv, a = a, b = b)
              },
              links = list(c("log", "identity", "sqrt"), "identity", "identity"),
              nam = "twlss",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = function(object, type = c("deviance", "pearson", "response")) {
                type <- match.arg(type)
                mu <- object$fitted.values[ , 1]
                p <- object$fitted.values[ , 2]
                ind <- p > 0
                ethi <- exp(-p[ind])
                ethni <- exp(p[!ind])
                p[ind] <- (b + a * ethi) / (1 + ethi)
                p[!ind] <- (b * ethni + a) / (ethni + 1)
                phi <- exp(object$fitted.values[ , 3])
                if (type == "pearson") {
                  rsd <- (object$y - mu) / sqrt(phi * mu^p)
                } else if (type == "response") {
                  rsd <- object$y - mu
                } else {
                  y1 <- object$y + (object$y == 0)
                  theta <- (y1^(1 - p) - mu^(1 - p)) / (1 - p)
                  kappa <- (object$y^(2 - p) - mu^(2 - p)) / (2 - p)
                  rsd <- sign(object$y - mu) * sqrt(pmax(2 * (object$y * theta - kappa) * object$prior.weights / phi, 0))
                }
                return( rsd )
              },
              postproc = mgcv::twlss()$postproc,
              # Here the arguments of twlss() are not used/important because twlss()$initialize is an
              # expression, not a closure. So the value of any variable it uses will depend on the
              # environment in which it is evaluated.
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled, weights){
                start <- NULL
                pen.reg <- penreg
                eval(twlss()$initialize)
                return( start )
              }
  )

  # Fixing the environment of all functions
  for(ii in 1:length(out)){
    if( class(out[[ii]]) == "function" ){
      environment(out[[ii]]) <- environment()
    }
  }

  return( out )
}
