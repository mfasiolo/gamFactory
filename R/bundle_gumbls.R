#'
#' Bundle for the Gumbel location-scale model
#'
#' @description Ported from [`mgcv::gumbls`]. The two linear predictors correspond,
#'              in order, to the mode \code{mu} (mandatory "identity" link) and
#'              \code{bt = log(beta)}, the log scale. The default link for \code{bt} is
#'              \code{"logea2(-7)"}, a numerically stabilized version of
#'              \code{"logea(exp(-7))"} (matching \code{mgcv::gumbls}'s own default of
#'              \code{b = -7}, see \link{make_link}); \code{"identity"} is also a valid,
#'              if less numerically robust, alternative.
#' @name bundle_gumbls
#' @rdname bundle_gumbls
#' @export
#'
bundle_gumbls <- function(){

  out <- list(np = 2,
              available_deriv = 4,
              llk = gamFactory::llk_gumbls,
              links = list("identity", "logea2(-7)"),
              nam = "gumbls",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = mgcv::gumbls()$residuals,
              rd = mgcv::gumbls()$rd,
              predict = mgcv::gumbls()$predict,
              postproc = mgcv::gumbls()$postproc,
              # Here the arguments of gumbls() are not used/important because gumbls()$initialize is an
              # expression, not a closure. So the value of any variable it uses will depend on the
              # environment in which it is evaluated.
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled, weights){
                start <- NULL
                pen.reg <- penreg
                eval(gumbls()$initialize)
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
