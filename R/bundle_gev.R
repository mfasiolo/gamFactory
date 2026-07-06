#'
#' Bundle for the Generalized Extreme Value (GEV) location-scale-shape model
#'
#' @description Ported from \code{mgcv::gevlss}. The three linear predictors correspond,
#'              in order, to the location \code{mu}, the log-scale \code{rho = log(sigma)}
#'              and the shape \code{xi} of the GEV distribution. The default link for
#'              \code{xi} is \code{"logitab(-1, 0.5)"}, bounding the shape in
#'              \code{(-1, 0.5)} (as in \code{mgcv::gevlss}'s own default logit-type
#'              link), which ensures a consistent estimator with finite variance.
#' @name bundle_gev
#' @rdname bundle_gev
#' @export
#'
bundle_gev <- function(){

  out <- list(np = 3,
              available_deriv = 4,
              llk = gamFactory::llk_gev,
              links = list(c("identity", "log"), "identity", "logitab(-1, 0.5)"),
              nam = "gev",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = mgcv::gevlss()$residuals,
              rd = mgcv::gevlss()$rd,
              postproc = mgcv::gevlss()$postproc,
              qf = function(p, mu, wt, scale, logp = FALSE) {
                mu <- as.matrix(mu)
                if (ncol(mu) == 1) { mu <- t(mu) }
                muE <- mu[ , 1, drop = TRUE]
                sigE <- exp(mu[ , 2, drop = TRUE])
                xiE <- mu[ , 3, drop = TRUE]
                if (logp) { p <- exp(p) }
                q <- muE + (sigE / xiE) * ((-log(p))^(-xiE) - 1)
                return( q )
              },
              cdf = function(q, mu, wt, scale, logp = FALSE) {
                mu <- as.matrix(mu)
                if (ncol(mu) == 1) { mu <- t(mu) }
                muE <- mu[ , 1, drop = TRUE]
                sigE <- exp(mu[ , 2, drop = TRUE])
                xiE <- mu[ , 3, drop = TRUE]
                z <- pmax(0, 1 + xiE * (q - muE) / sigE)
                p <- -z^(-1/xiE)
                if (!logp) { p <- exp(p) }
                return( p )
              },
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled, weights){
                start <- NULL
                pen.reg <- penreg
                # Here the arguments of gevlss() are not used/important because gevlss()$initialize is an expression, 
                # not a closure. So the value of any variable it uses will depend on the environment in which 
                # it is evaluated. 
                eval(gevlss()$initialize)
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
