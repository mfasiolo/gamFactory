#'
#' Ingredients for two states Markov model
#' 
#' @description XXX.
#' @name bundleMARKOV
#' @rdname bundleMARKOV
#' @export bundleMARKOV
#' @examples
#' 
#' library(gamFactory)
#' 
#' # Simulate some data from Markov chain with two states (0, 1)
#' n <- 5000
#' x = sort(runif(n, -1, 1))
#' 
#' # Probabilities of staying in state 0 and 1
#' p1 <- exp(x^2) / (1 + exp(x^2))
#' p2 <- exp(-x - 0.6*x^2) / (1 + exp(-x - 0.6*x^2))
#' 
#' # Simulate observations 
#' y0 <- numeric(n)
#' y0[1] <- 1
#' for(ii in 2:n){
#'   if(y0[ii-1] == 0) y0[ii] <- rbinom(1, 1, 1 - p1[ii])
#'   if(y0[ii-1] == 1) y0[ii] <- rbinom(1, 1, p2[ii])
#' }
#' 
#' myDat <- data.frame(x = x, y = y0)
#' 
#' # Plot data and probabilities
#' par(mfrow = c(2, 2))
#' plot(myDat)
#' plot(x, p1)
#' plot(x, p2)
#' 
#' # Create Markov model and fit it
#' myMRK <- buildGAMLSS( bundleMARKOV )
#' fit <- gam(list(y ~ s(x), ~ s(x)), data = myDat, family = myMRK)
#' 
#' # Compare true and estimated probabilities
#' plot(fit, select = 1, 
#'      trans = function(.x) exp(coef(fit)[1] + .x) / (1 + exp(coef(fit)[1] + .x))) 
#' lines(x, p1, col = 2)
#' 
#' plot(fit, select = 2, 
#'      trans = function(.x) exp(coef(fit)[11] + .x) / (1 + exp(coef(fit)[11] + .x)))
#' lines(x, p2, col = 2)
#' 
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

