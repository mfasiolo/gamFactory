#'
#' The gaussian family
#' 
#' @name fam_gaussian
#' @rdname fam_gaussian
#' @export
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
fam_gaussian <- function(){
  
  bundle <- bundle_gaussian()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}

