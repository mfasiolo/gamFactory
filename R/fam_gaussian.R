#'
#' The Gaussian location-scale family
#'
#' @description \code{gam}/\code{gam_nl} should be called with a list of 2 formulae: the
#'              first specifies the response and the linear predictor for the mean
#'              \code{mu}; the second (one sided) specifies the linear predictor for
#'              \code{tau = 1 / sigma}, the reciprocal of the standard deviation.
#' @name fam_gaussian
#' @rdname fam_gaussian
#' @export fam_gaussian
#' @examples
#' library(gamFactory)
#' 
#'  
#' library(MASS)
#' 
#' # Create gaussian family and fit it
#' fit <- gam_nl(list(accel ~ s(times,k=20,bs="ad"),
#'                    ~ s(times)), data = mcycle, family = fam_gaussian())
#' plot(fit, pages = 1, scale = FALSE)
#' 
#' # Compare true (black) and estimated (red) quantiles
#' qu <- 0.9
#' plot(mcycle, col = "grey")
#' lines(mcycle$times, fit$family$qf(qu, mu = fit$fitted.values), col = 4)
#' 
#' fit_2 <- gam(list(accel ~ s(times,k=20,bs="ad"),
#'                   ~ s(times)), data = mcycle, family = fam_gaussian(),
#'              optimizer = c("outer", "bfgs"))
#' 
#' err <- max(abs(fit$fitted.values - fit_2$fitted.values))
#' if(err > 1e-6){
#'   stop("Discrepancy between gam and gam_nl")
#' }
fam_gaussian <- function(link = NULL){
  
  bundle <- bundle_gaussian()
  
  fam <- build_family(bundle, link = link)()
  
  return(fam)
  
}

