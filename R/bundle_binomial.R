#'
#' Bundle for Gaussian regression model
#' 
#' @name bundle_binomial
#' @rdname bundle_binomial
#' @export
#'
bundle_binomial <- function(n){
  force(n)
  # We want to fix n (the "size" of the binomial) at this stage
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
                fam_nam <- "binomial"
                y <- drop(object$y)
                if(length(n) != length(y)){
                  if(length(n) == 1){
                    n <- rep(n, length(y))
                  } else {
                    stop("length(n) != length(y)")
                  }
                }
                # We need to do something special here to deal with n
                poi_fam <- fam <- do.call(fam_nam, list())
                fam$dev.resid <- function(y, mu, wt){
                 wt * poi_fam$dev.resid(y/n, mu, n)
                }
                fam$var <- function(mu){
                 n * poi_fam$var(mu) 
                }
                fam$exp_val <- function(mu){
                  n * mu
                }
                # Make sure these function use "n" from the current environment
                environment(fam$dev.resid) <- environment(fam$var) <- 
                  environment(fam$exp_val) <- environment()
                r <- .resid_exp_fam(object = object, type = type, fam = fam)
                return( r )
              },
              rd = function(mu, wt, scale) {
                return( rbinom(nrow(mu), n, mu) )
              },
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
  
  # Fixing the environment of all functions
  for(ii in 1:length(out)){
    if( class(out[[ii]]) == "function" ){
      environment(out[[ii]]) <- environment()
    }
  }
  
  return( out )
}

