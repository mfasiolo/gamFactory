#'
#' Ingredients for Generalized Pareto family
#' 
#' @description XXX.
#' @name bundle_gpd
#' @rdname bundle_gpd
#' @export bundle_gpd
bundle_gpd <- function(){
  out <- list(np = 2,
              available_deriv = 3,
              llk = gamFactory::llk_gpd,
              links = list(c("log", "sqrt"), "logitab(0.001, 0.5)"), 
              nam = "gpd", 
              bundle_nam = as.character(match.call()[[1]]),
              residuals = function(object, type = c("deviance", "pearson", "response")) {
                type <- match.arg(type)
                xi <- object$fitted[ , 2] 
                sig <- object$fitted[ , 1] / (1 + xi) # sig = phi / (1 + xi)
                rsd <- object$y - sig / (1 - xi)      # E(y) = sig / (1 - xi)
                if (type=="response"){ 
                  return( rsd )
                } else {
                  return( rsd / (sig / ((1-xi)*sqrt(1-2*xi))) ) # sd(y) = sig / [ (1 - xi) * sqrt(1 - 2 * xi)   
                }
              }, 
              qf = function(p, mu, wt, scale) {
                mu <- as.matrix(mu)
                if( ncol(mu) == 1 ){ mu <- t(mu) }
                xi <- mu[ , 2, drop = TRUE]
                sig <- mu[ , 1, drop = TRUE] / (1 + xi)
                q <- sig * ((1-p)^(-xi) - 1) / xi
                return( q )
              }, 
              cdf = function(q, mu, wt, scale, logp) {
                mu <- as.matrix(mu)
                if( ncol(mu) == 1 ){ mu <- t(mu) }
                xi <- mu[ , 2, drop = TRUE]
                sig <- mu[ , 1, drop = TRUE] / (1 + xi)
                p <- log1p( - (1 + xi * q / sig)^(-1/xi) ) 
                if( !logp ){ p <- exp(p)  }
                return( p )
              }, 
              rd = function(mu, wt, scale) { 
                mu <- as.matrix(mu)
                if( ncol(mu) == 1 ){ mu <- t(mu) }
                xi <- mu[ , 2, drop = TRUE]
                sig <- mu[ , 1, drop = TRUE] / (1 + xi)
                
                n <- length(xi)
                r <- sig * ((1-runif(n))^(-xi) - 1) / xi
                
                return( r )
              }, 
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                ## start out with xi close to zero. If xi==0 then
                ## mean is sigma = phi. Idea is to regress g(y) on model matrix for mean.
                ## Note that appropriate E scaling
                ## for full calculation may be inappropriate for initialization 
                ## which is basically penalizing something different here.
                ## best we can do here is to use E only as a regularizer.
                
                n <- rep(1, nobs)
                
                ## should E be used unscaled or not?..
                use.unscaled <- if (!is.null(attr(E, "use.unscaled"))){ TRUE } else { FALSE }
                
                lpi <- attr(x,"lpi")
                start <- rep(0, ncol(x))
                
                #### Get scale phi = sigma = E(y) (if xi == 0)
                yt1 <- if (family$link[[1]]=="identity"){ 
                  y 
                } else {
                  family$linfo[[1]]$linkfun(abs(y) + max(y)*1e-7)
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
                
                #### Get shape xi to be the equivalent of nearly 0
                x1 <-  x[ , lpi[[2]], drop=FALSE]
                startji <- qr.coef(qr(x1), c(rep(family$linfo[[2]]$linkfun(1e-3),nrow(x1))))   
                startji[!is.finite(startji)] <- 0
                start[lpi[[2]]] <- startji
                
                return( start )
              }
  )
  
  return(out)
}