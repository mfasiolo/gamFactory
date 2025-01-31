#'
#' Bundle for Sinh-Arsinh (shash) regression model
#' 
#' @name bundle_shash
#' @rdname bundle_shash
#' @export
#'
bundle_shash <- function(){
  out <- list(np = 4,
              available_deriv = 3,
              llk = gamFactory:::llk_shash,
              links = list(c("identity", "inverse", "log", "sqrt"), "logea(0.01)", "identity", "identity"), 
              nam = "shash",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = mgcv::shash()$residuals,
              rd = mgcv::shash()$rd,
              cdf = mgcv::shash()$cdf,
              qf = function(p, mu, wt, scale, logp = FALSE) {
                mu <- as.matrix(mu)
                if(ncol(mu)==1){ mu <- t(mu) }
                muE <- mu[ , 1, drop = TRUE]
                sigE <- exp(mu[ , 2, drop = TRUE])
                epsE <- mu[ , 3, drop = TRUE]
                delE <- exp(mu[ , 4, drop = TRUE])
                q <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(p, log.p = logp)) + (epsE/delE))
                return(q)
              },
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                start <- NULL
                pen.reg <- penreg
                eval(shash()$initialize)
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

