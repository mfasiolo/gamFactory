#'
#' Bundle for the Zero-Inflated Poisson (ZIP) location-scale model
#'
#' @description Ported from [`mgcv::ziplss`]. The two linear predictors correspond,
#'              in order, to \code{g = log(lambda)} (the Poisson rate of the "at risk"
#'              component) and \code{eta} (the zero-inflation probability, related via
#'              \code{p = 1 - exp(-exp(eta))}). Both linear predictors only support the
#'              "identity" link, as in \code{mgcv::ziplss}.
#' @name bundle_ziplss
#' @rdname bundle_ziplss
#' @export
#'
bundle_ziplss <- function(){

  out <- list(np = 2,
              available_deriv = 4,
              llk = gamFactory::llk_ziplss,
              links = list("identity", "identity"),
              nam = "ziplss",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = mgcv::ziplss()$residuals,
              rd = mgcv::ziplss()$rd,
              predict = mgcv::ziplss()$predict,
              postproc = mgcv::ziplss()$postproc,
              # Here the arguments of ziplss() are not used/important because ziplss()$initialize is an
              # expression, not a closure. So the value of any variable it uses will depend on the
              # environment in which it is evaluated.
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled, weights){
                start <- NULL
                pen.reg <- penreg
                eval(ziplss()$initialize)
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
