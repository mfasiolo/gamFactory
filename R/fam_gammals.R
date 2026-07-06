#'
#' The Gamma location-scale family
#'
#' @description Ported from [mgcv::gammals]. \code{gam} should be called with a list
#'              of 2 formulae: the first specifies the response and the linear predictor
#'              for \code{mu = log(mean)}; the second (one sided) specifies the linear
#'              predictor for \code{th = log(phi)}, the log dispersion
#'              (\code{phi = 1 / shape}).
#' @name fam_gammals
#' @rdname fam_gammals
#' @export fam_gammals
#' @examples
#' library(gamFactory)
#'
#' ## Simulate data from a Gamma location-scale model
#' set.seed(1)
#' n <- 1000
#' f1 <- function(x) 1 + sin(2 * pi * x)
#' f2 <- function(x) 1 - x
#'
#' x1 <- runif(n); x2 <- runif(n)
#' mu <- f1(x1)
#' phi <- exp(f2(x2))
#'
#' y <- rgamma(n, shape = 1 / phi, scale = mu * phi)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Fit model
#' fit <- gam(list(y ~ s(x1), ~ s(x2)), data = dat, family = fam_gammals())
#' plot(fit, pages = 1, scale = 0)
#'
#' muE <- exp(fit$fitted[ , 1])
#' thE <- fit$fitted[ , 2]
#'
#' par(mfrow = c(1, 2))
#' plot(x1, muE); lines(sort(x1), mu[order(x1)], col = 2)
#' plot(x2, thE); lines(sort(x2), log(phi)[order(x2)], col = 2)
#' 
#' # Compare with mgcv
#' fit_mgcv <- gam(list(y ~ s(x1), ~ s(x2)), data = dat, family = gammals())
#' 
#' par(mfrow = c(1, 2))
#' plot(exp(fit_mgcv$fitted.value[ , 1]), muE); abline(0, 1, col = 2)
#' plot(fit_mgcv$fitted.value[ , 2], thE); abline(0, 1, col = 2)
#'
fam_gammals <- function(link = NULL){

  bundle <- bundle_gammals()

  fam <- build_family(bundle, link = link)()

  return(fam)

}
