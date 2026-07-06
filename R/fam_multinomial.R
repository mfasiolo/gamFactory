#'
#' The multinomial family
#'
#' @description Ported from [mgcv::multinom], see the corresponding documentation for details.
#' @param K number of linear predictors, i.e. number of categories minus one.
#' @name fam_multinomial
#' @rdname fam_multinomial
#' @export fam_multinomial
#' @examples
#' library(gamFactory)
#'
#' set.seed(0)
#' n <- 500; K <- 2
#' x <- runif(n)
#' eta <- cbind(2 * x - 1, 1 - 2 * x)
#' p <- cbind(1, exp(eta)); p <- p / rowSums(p)
#' y <- apply(p, 1, function(pr) sample(0:K, 1, prob = pr))
#' dat <- data.frame(y = y, x = x)
#'
#' fit <- gam(list(y ~ s(x), ~ s(x)), data = dat, family = fam_multinomial(K = K))
#' plot(fit, pages = 1)
#'
fam_multinomial <- function(K, link = NULL){

  bundle <- bundle_multinomial(K)

  fam <- build_family(bundle, link = link)()

  return(fam)

}
