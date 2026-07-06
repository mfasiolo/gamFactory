#'
#' Bundle for the Gamma location-scale model
#'
#' @description Ported from [mgcv::gammals]. The two linear predictors correspond,
#'              in order, to \code{mu = log(mean)} (mandatory "identity" link, but the
#'              value it maps to is internally the log of the mean, as in
#'              \code{mgcv::gammals}) and \code{th = log(phi)}, the log dispersion
#'              (\code{phi = 1 / shape}). The default link for \code{th} is
#'              \code{"logea2(-7)"}, a numerically stabilized version of
#'              \code{"logea(exp(-7))"} (matching \code{mgcv::gammals}'s own default of
#'              \code{b = -7}, see \link{make_link}); \code{"identity"} is also a valid,
#'              if less numerically robust, alternative.
#' @name bundle_gammals
#' @rdname bundle_gammals
#' @export
#'
bundle_gammals <- function(){

  out <- list(np = 2,
              available_deriv = 4,
              llk = gamFactory::llk_gammals,
              links = list("identity", c("logea2(-7)", "identity")),
              nam = "gammals",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = mgcv::gammals()$residuals,
              rd = mgcv::gammals()$rd,
              predict = mgcv::gammals()$predict,
              postproc = mgcv::gammals()$postproc,
              # Here the arguments of gammals() are not used/important because gammals()$initialize is an
              # expression, not a closure. So the value of any variable it uses will depend on the
              # environment in which it is evaluated.
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled, weights){
                start <- NULL
                pen.reg <- penreg
                eval(gammals()$initialize)
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
