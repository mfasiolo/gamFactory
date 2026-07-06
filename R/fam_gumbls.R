#'
#' The Gumbel location-scale family
#'
#' @description Ported from [`mgcv::gumbls`]. \code{gam} should be called with a list
#'              of 2 formulae: the first specifies the response and the linear predictor
#'              for the mode \code{mu}; the second (one sided) specifies the linear
#'              predictor for \code{bt = log(beta)}, the log scale.
#' @name fam_gumbls
#' @rdname fam_gumbls
#' @export fam_gumbls
#' @examples
#' library(gamFactory)
#'
#' ## Simulate data from a Gumbel location-scale model
#' set.seed(1)
#' n <- 1000
#' f1 <- function(x) 2 + 2 * sin(2 * pi * x)
#' f2 <- function(x) 1 - x
#'
#' x1 <- runif(n); x2 <- runif(n)
#' mu <- f1(x1)
#' beta <- exp(f2(x2))
#'
#' u <- runif(n)
#' y <- mu - beta * log(-log(u))
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Fit model
#' fit <- gam(list(y ~ s(x1), ~ s(x2)), data = dat, family = fam_gumbls())
#' plot(fit, pages = 1, scale = 0)
#'
#' muE <- fit$fitted[ , 1]
#' btE <- fit$fitted[ , 2]
#'
#' par(mfrow = c(1, 2))
#' plot(x1, muE); lines(sort(x1), mu[order(x1)], col = 2)
#' plot(x2, btE); lines(sort(x2), log(beta)[order(x2)], col = 2)
#'
fam_gumbls <- function(link = NULL){

  bundle <- bundle_gumbls()

  fam <- build_family(bundle, link = link)()

  return(fam)

}
