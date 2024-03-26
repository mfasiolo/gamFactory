#'
#' Bundle for Gaussian regression model
#' 
#' @name bundle_binomial
#' @rdname bundle_binomial
#' @export
#'
bundle_binomial <- function(n){
  force(n)
  .llk_wrap <- function(...){
    llk_binomial(n = n, ...)
  }
  out <- list(np = 1,
              available_deriv = 4,
              llk = .llk_wrap,
              links = list(c("logit", "probit")), 
              nam = "binomial",
              store = list("n" = n),
              bundle_nam = as.character(match.call()[[1]]),
              residuals = function(object, type=c("deviance", "pearson", "response")) {
                type <- match.arg(type)
                rsd <- object$y-object$fitted[,1]
                #    if (type=="response") return(rsd) else
                #      return((rsd*object$fitted[,2])) ## (y-mu)/sigma 
              },
              # rd = function(mu, wt, scale) {
              #   return( rnorm(nrow(mu), mu[ , 1], sqrt(scale/wt)/mu[ , 2]) )
              # },
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                
                n <- family$store$n
                
                ## should E be used unscaled or not?..
                use.unscaled <- if (!is.null(attr(E, "use.unscaled"))){ TRUE } else { FALSE }
                
                lpi <- attr(x,"lpi")
                start <- rep(0, ncol(x))
                
                yt1 <- if (family$link[[1]]=="identity"){ 
                  y 
                } else {
                  family$linfo[[1]]$linkfun(pmax(pmin(y/n, 1-1e-3), 1e-3))
                }
                x1 <- x[ , lpi[[1]], drop=FALSE]
                e1 <- E[ , lpi[[1]], drop=FALSE] ## square root of total penalty
                
                if (use.unscaled) {
                  qrx <- qr( rbind(x1, e1) )
                  x1 <- rbind(x1, e1)
                  startji <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
                  startji[!is.finite(startji)] <- 0       
                } else {
                  startji <- penreg(x1, e1, yt1)
                }
                start[lpi[[1]]] <- startji
                
                return( start )
              }
  )
  return( out )
}

