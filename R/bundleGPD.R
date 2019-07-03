#'
#' Ingredients for Generalized Pareto family
#' 
#' @description XXX.
#' @name bundleGPD
#' @rdname bundleGPD
#' @export bundleGPD
#' @examples 
#' library(gamFactory)
#'
#' # Simulate some data
#' n <- 2000
#' myDat <- data.frame(x = runif(n, -1, 1))
#' myDat$y <-  rexp(n, 1 / (0.3 + abs(myDat$x)^2))
#' plot(myDat)
#' 
#' # Create GPD family and fit it
#' library(mgcViz)
#' myGPD <- buildGAMLSS( bundleGPD )
#' fit <- gamV(list(y ~ s(x), ~ s(x)), data = myDat, family = myGPD)
#' print(plot(fit), pages = 1)
#' 
#' # Compare true (black) and estimated (red) quantiles 
#' qu <- 0.95
#' xseq <- seq(-1, 1, length.out = 1000)
#' plot(myDat, col = "grey")
#' lines(xseq, qexp(qu, 1 / (0.3 + abs(xseq)^2)), col = 2)
#' lines(xseq, fit$family$qf(qu, 
#'                           mu = predict(fit, newdata = data.frame(x = xseq), 
#'                                        type = "response")), 
#'       col = 3)
#'
bundleGPD <- list(np = 2,
                  llk = function(y, param, deriv) logLikGPD(y)$derObj(param, deriv),
                  links = list(c("log", "sqrt"), "logitab(-0.5, 0.5)"), 
                  nam = "gpd", 
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
                  initialize = expression({
                    ## start out with xi close to zero. If xi==0 then
                    ## mean is sigma = phi. Idea is to regress g(y) on model matrix for mean.
                    ## Note that appropriate E scaling
                    ## for full calculation may be inappropriate for initialization 
                    ## which is basically penalizing something different here.
                    ## best we can do here is to use E only as a regularizer.
                    n <- rep(1, nobs)
                    ## should E be used unscaled or not?..
                    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))){ TRUE } else { FALSE }
                    if ( is.null(start) ) {
                      jj <- attr(x,"lpi")
                      start <- rep(0, ncol(x))
                      
                      #### Get scale phi = sigma = E(y) (if xi == 0)
                      yt1 <- if (family$link[[1]]=="identity"){ 
                        y 
                      } else {
                        family$linfo[[1]]$linkfun(abs(y) + max(y)*1e-7)
                      }
                      x1 <- x[ , jj[[1]], drop=FALSE]
                      e1 <- E[ , jj[[1]], drop=FALSE] ## square root of total penalty
                      
                      if (use.unscaled) {
                        qrx <- qr( rbind(x1, e1) )
                        x1 <- rbind(x1, e1)
                        startji <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
                        startji[!is.finite(startji)] <- 0       
                      } else {
                        startji <- pen.reg(x1, e1, yt1)
                      }
                      start[jj[[1]]] <- startji
                      
                      #### Get shape xi to be the equivalent of nearly 0
                      x1 <-  x[ , jj[[2]], drop=FALSE]
                      startji <- qr.coef(qr(x1), c(rep(family$linfo[[2]]$linkfun(1e-3),nrow(x1))))   
                      startji[!is.finite(startji)] <- 0
                      start[jj[[2]]] <- startji
                    }
                  }) 
)