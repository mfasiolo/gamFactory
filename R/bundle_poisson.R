#'
#' Bundle for Gaussian regression model
#' 
#' @name bundle_poisson
#' @rdname bundle_poisson
#' @export bundle_poisson
#'
bundle_poisson <- function(){
  out <- list(np = 1,
              available_deriv = 4,
              llk = gamFactory::llk_poisson,
              links = list(c("log", "sqrt")), 
              nam = "poisson",
              bundle_nam = as.character(match.call()[[1]]),
              residuals = function(object, type=c("deviance", "pearson", "response")) {
                type <- match.arg(type)
                fam <- do.call("poisson", list())
                r <- .resid_exp_fam(object = object, type = type, fam = fam)
                return( r )
              },
              rd = function(mu, wt, scale) {
                return( rpois(nrow(mu), mu) )
              },
              initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                
                n <- rep(1, nobs)
                
                ## should E be used unscaled or not?..
                use.unscaled <- if (!is.null(attr(E, "use.unscaled"))){ TRUE } else { FALSE }
                
                lpi <- attr(x,"lpi")
                start <- rep(0, ncol(x))
                
                #### Get scale phi = sigma = E(y) (if xi == 0)
                yt1 <- if (family$link[[1]]=="identity"){ 
                  y 
                } else {
                  family$linfo[[1]]$linkfun(as.double(y) + 1e-3)
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
  
  # Fixing the environment of all functions
  for(ii in 1:length(out)){
    if( class(out[[ii]]) == "function" ){
      environment(out[[ii]]) <- environment()
    }
  }
  
  return( out )
}

