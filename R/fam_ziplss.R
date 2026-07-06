#'
#' The Zero-Inflated Poisson (ZIP) location-scale family
#'
#' @description Ported from [`mgcv::ziplss`]. \code{gam} should be called with a list
#'              of 2 formulae: the first specifies the response and the linear predictor
#'              for \code{g = log(lambda)}, the Poisson rate of the "at risk" component;
#'              the second (one sided) specifies the linear predictor for the
#'              zero-inflation probability, related via \code{p = 1 - exp(-exp(eta))}.
#' @name fam_ziplss
#' @rdname fam_ziplss
#' @export fam_ziplss
#' @examples
#' library(gamFactory)
#'
#' ## Simulate data from a ZIP location-scale model
#' set.seed(1)
#' n <- 1000
#' f1 <- function(x) 2 * sin(pi * x)
#' f2 <- function(x) 1 + x - 2 * x^2
#'
#' x1 <- runif(n); x2 <- runif(n)
#' g <- 1 + f1(x1)     # linear predictor for log(lambda)
#' eta <- f2(x2)        # linear predictor for the zero-inflation probability
#'
#' lambda <- exp(g)
#' p <- 1 - exp(-exp(eta))
#' y <- ifelse(p > runif(n), rpois(n, lambda), 0)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Fit model
#' fit <- gam(list(y ~ s(x1), ~ s(x2)), data = dat, family = fam_ziplss())
#' plot(fit, pages = 1, scale = 0)
#'
#' gE <- fit$fitted[ , 1]
#' etaE <- fit$fitted[ , 2]
#'
#' par(mfrow = c(1, 2))
#' plot(x1, gE); lines(sort(x1), g[order(x1)], col = 2)
#' plot(x2, etaE); lines(sort(x2), eta[order(x2)], col = 2)
#'
fam_ziplss <- function(link = NULL){

  bundle <- bundle_ziplss()

  fam <- build_family(bundle, link = link)()

  return(fam)

}
