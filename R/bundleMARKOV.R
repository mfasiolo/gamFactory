#'
#' Ingredients for two states Markov model
#' 
#' @description XXX.
#' @name bundleMARKOV
#' @rdname bundleMARKOV
#' @export bundleMARKOV
#' @examples 
#' library(gamFactory)
#' 
#' # Simulate some data
#' n <- 2000
#' myDat <- data.frame(x = sort(runif(n, -1, 1)))
#' p1 <- pnorm(myDat$x, 0, 0.7)[-1]
#' p2 <- pnorm(-myDat$x, 0, 0.7)[-1]
#' 
#' y0 <- numeric(n)
#' y0[1] <- 1
#' for(ii in 2:n){
#'   if(y0[ii-1] == 0) y0[ii] <- rbinom(1, 1, 1 - p1[ii])
#'   if(y0[ii-1] == 1) y0[ii] <- rbinom(1, 1, 1 - p1[ii])
#' }
#' 
#' myDat$y <-  y0
#' plot(myDat)
#' 
#' # Create GPD family and fit it
#' library(mgcViz)
#' myMRK <- buildGAMLSS( bundleMARKOV )
#' fit <- gamV(list(y ~ s(x), ~ s(x)), data = myDat, family = myMRK)
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
#' #' 
bundleMARKOV <- list(np = 2,
                  availableDeriv = 3,
                  llk = function(y, param, deriv, ...) logLikMARKOV(y)$derObj(param, deriv),
                  links = list(c("logit"), c("logit")), 
                  nam = "markov", 
                  residuals = function(object, type = c("deviance", "pearson", "response")) {
                    type <- match.arg(type)
                    # xi <- object$fitted[ , 2] 
                    # sig <- object$fitted[ , 1] / (1 + xi) # sig = phi / (1 + xi)
                    # rsd <- object$y - sig / (1 - xi)      # E(y) = sig / (1 - xi)
                    if (type=="response"){ 
                      return( object$y )
                    } else {
                      return( object$y ) # sd(y) = sig / [ (1 - xi) * sqrt(1 - 2 * xi)   
                    }
                  }, 
                  initFun = function(y, nobs, E, x, family){
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

                    return( start )
                  }
)

