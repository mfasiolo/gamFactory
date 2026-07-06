#'
#' The Generalized Extreme Value (GEV) location-scale-shape family
#'
#' @description Ported from [`mgcv::gevlss`]. \code{gam} should be called with a list
#'              of 3 formulae: the first specifies the response and the linear predictor
#'              for the location parameter, \code{mu}; the second (one sided) specifies
#'              the linear predictor for the log scale parameter, \code{log(sigma)}; the
#'              third (one sided) specifies the linear predictor for the shape parameter,
#'              \code{xi}.
#' @name fam_gev
#' @rdname fam_gev
#' @export fam_gev
#' @examples
#' library(gamFactory)
#'
#' ## Simulate data from a GEV location-scale-shape model
#' set.seed(1)
#' n <- 1000
#' f0 <- function(x) 2 * sin(pi * x)
#' f1 <- function(x) exp(2 * x)
#' f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
#'
#' x0 <- runif(n); x1 <- runif(n); x2 <- runif(n)
#' mu <- f2(x2)
#' rho <- f0(x0)
#' xi <- (f1(x1) - 4) / 9
#'
#' # Inverse cdf of the GEV distribution
#' Fi.gev <- function(z, mu, sigma, xi) {
#'   xi[abs(xi) < 1e-8] <- 1e-8
#'   mu + ((-log(z))^(-xi) - 1) * sigma / xi
#' }
#' y <- Fi.gev(runif(n), mu, exp(rho), xi)
#' dat <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2)
#'
#' ## Fit model
#' fit <- gam(list(y ~ s(x2), ~ s(x0), ~ s(x1)), data = dat, family = fam_gev())
#' plot(fit, pages = 1, scale = 0)
#'
#' muE <- fit$fitted[ , 1]
#' rhoE <- fit$fitted[ , 2]
#' xiE <- fit$fitted[ , 3]
#'
#' par(mfrow = c(2, 2))
#' plot(x2, muE); lines(sort(x2), mu[order(x2)], col = 2)
#' plot(x0, rhoE); lines(sort(x0), rho[order(x0)], col = 2)
#' plot(x1, xiE); lines(sort(x1), xi[order(x1)], col = 2)
#'
fam_gev <- function(link = NULL){

  bundle <- bundle_gev()

  fam <- build_family(bundle, link = link)()

  return(fam)

}
