#'
#' The Tweedie location-scale-shape family
#'
#' @description Ported from [`mgcv::twlss`]. \code{gam} should be called with a list
#'              of 3 formulae: the first specifies the response and the linear predictor
#'              for the mean \code{mu}; the second (one sided) specifies the linear
#'              predictor for \code{theta}, which determines the Tweedie power parameter
#'              via \code{p = (a + b * exp(theta)) / (1 + exp(theta))}; the third (one
#'              sided) specifies the linear predictor for \code{rho = log(phi)}, the log
#'              dispersion.
#' @param a lower limit parameter (> 1) used in the definition of the Tweedie power
#'          parameter \code{p} from \code{theta}. See \code{mgcv::ldTweedie}.
#' @param b upper limit parameter (< 2) used in the definition of the Tweedie power
#'          parameter \code{p} from \code{theta}. See \code{mgcv::ldTweedie}.
#' @name fam_twlss
#' @rdname fam_twlss
#' @export fam_twlss
#' @examples
#' library(gamFactory)
#'
#' ## Simulate data from a Tweedie location-scale-shape model
#' set.seed(1)
#' n <- 1000
#' f1 <- function(x) 2 + 2 * sin(pi * x)
#' f2 <- function(x) 1 - x
#'
#' x1 <- runif(n); x2 <- runif(n)
#' mu <- exp(f1(x1) / 3)
#' phi <- exp(f2(x2))
#' p <- 1.5
#'
#' y <- rTweedie(mu, p = p, phi = phi)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Fit model (theta held constant, i.e. p is not estimated by a smooth)
#' fit <- gam(list(y ~ s(x1), ~ 1, ~ s(x2)), data = dat, family = fam_twlss())
#' plot(fit, pages = 1, scale = 0)
#'
#' muE <- fit$fitted[ , 1]
#' rhoE <- fit$fitted[ , 3]
#'
#' par(mfrow = c(1, 2))
#' plot(x1, muE); lines(sort(x1), mu[order(x1)], col = 2)
#' plot(x2, rhoE); lines(sort(x2), log(phi)[order(x2)], col = 2)
#'
fam_twlss <- function(link = NULL, a = 1.01, b = 1.99){

  bundle <- bundle_twlss(a = a, b = b)

  fam <- build_family(bundle, link = link)()

  return(fam)

}
