#'
#' Bundle for Gaussian regression model
#' 
#' @name bundle_gaussian
#' @rdname bundle_gaussian
#' @export
#'
bundle_gaussian <- function(){
  out <- list(np = 2,
              available_deriv = 4,
              llk = gamFactory::llk_gaussian,
              links = list(c("identity", "inverse", "log", "sqrt"), "loga(0.01)"), 
              nam = "gaussian",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = function(object, type=c("deviance", "pearson", "response")) {
                type <- match.arg(type)
                rsd <- object$y-object$fitted[,1]
                if (type=="response") return(rsd) else
                  return((rsd*object$fitted[,2])) ## (y-mu)/sigma 
              },
              rd = function(mu, wt, scale) {
                return( rnorm(nrow(mu), mu[ , 1], sqrt(scale/wt)/mu[ , 2]) )
              },
              cdf = function(q, mu, wt, scale, logp = FALSE) {
                return( pnorm(q, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp) )
              },
              qf = function(p, mu, logp = FALSE) {
                return( qnorm(p, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp) )
              },
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                start <- NULL
                pen.reg <- penreg
                eval(gaulss()$initialize)
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


# Regress g(y) on model matrix for mean and corresponding log absolute residuals on 
# the model matrix for log(sigma) 
# n <- rep(1, nobs)
# use.unscaled <- if (!is.null(unscaled)) TRUE else FALSE
# if (!is.null(offset)) offset[[3]] <- 0
# yt1 <- if (family$link[[1]]=="identity") y else
#   family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
# if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
# start <- rep(0, ncol(x))
# x1 <- x[,jj[[1]],drop=FALSE]
# e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
# if (use.unscaled) {
#   qrx <- qr(rbind(x1,e1))
#   x1 <- rbind(x1,e1)
#   startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
#   startji[!is.finite(startji)] <- 0
# } else {
#   startji <- penreg(x1,e1,yt1)
# }
# start[jj[[1]]] <- startji
# lres1 <- log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])))
# if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
# x1 <- x[,jj[[2]],drop=FALSE]
# e1 <- E[,jj[[2]],drop=FALSE]
# if (use.unscaled) {
#   x1 <- rbind(x1,e1)
#   startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))
#   startji[!is.finite(startji)] <- 0
# } else {
#   startji <- penreg(x1,e1,lres1)
# }
# start[jj[[2]]] <- startji
# return( start )