#'
#' The Sinh-Arsinh (shash) regression model
#' 
#' @name fam_shash
#' @rdname fam_shash
#' @export fam_shash
#' @examples 
#' ###############
#' # [1] Example with standard effects 
#' ###############
#' ##  Simulate some data from shash
#' set.seed(847)
#' n <- 1000
#' x <- seq(-4, 4, length.out = n)
#' 
#' X <- cbind(1, x, x^2)
#' beta <- c(4, 1, 1)
#' mu <- X %*% beta 
#' 
#' sigma =  .5+0.4*(x+4)*.5            # Scale
#' eps = 2*sin(x)                      # Skewness
#' del = 1 + 0.2*cos(3*x)              # Kurtosis
#' 
#' dat <-  mu + (del*sigma)*sinh((1/del)*asinh(qnorm(runif(n))) + (eps/del))
#' dataf <- data.frame(cbind(dat, x))
#' names(dataf) <- c("y", "x")
#' plot(x, dat, xlab = "x", ylab = "y")
#' 
#' my_fam <- build_family(bundle_shash())
#' 
#' ## Fit model
#' fit <- gam(list(y ~ s(x), # <- model for location 
#'                 ~ s(x),   # <- model for log-scale
#'                 ~ s(x),   # <- model for skewness
#'                 ~ s(x, k = 20)), # <- model for log-kurtosis
#'            data = dataf, 
#'            family = my_fam, # <- new family
#'            )
#' 
#' muE <- fit$fitted[ , 1]
#' sigE <- exp(fit$fitted[ , 2])
#' epsE <- fit$fitted[ , 3]
#' delE <- exp(fit$fitted[ , 4])
#' 
#' par(mfrow = c(2, 2))
#' plot(x, muE, type = 'l', ylab = expression(mu(x)), lwd = 2)
#' lines(x, mu, col = 2, lty = 2, lwd = 2)
#' legend("top", c("estimated", "truth"), col = 1:2, lty = 1:2, lwd = 2)
#' 
#' plot(x, sigE, type = 'l', ylab = expression(sigma(x)), lwd = 2)
#' lines(x, sigma, col = 2, lty = 2, lwd = 2)
#' 
#' plot(x, epsE, type = 'l', ylab = expression(epsilon(x)), lwd = 2)
#' lines(x, eps, col = 2, lty = 2, lwd = 2)
#' 
#' plot(x, delE, type = 'l', ylab = expression(delta(x)), lwd = 2)
#' lines(x, del, col = 2, lty = 2, lwd = 2)
#' 
#' ## Plotting true and estimated conditional density
#' par(mfrow = c(1, 1))
#' plot(x, dat, pch = '.', col = "grey", ylab = "y", ylim = c(-35, 70))
#' for(qq in c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)){
#'   est <- fit$family$qf(p=qq, mu = fit$fitted)
#'   true <- mu + (del * sigma) * sinh((1/del) * asinh(qnorm(qq)) + (eps/del))
#'   lines(x, est, type = 'l', col = 1, lwd = 2)
#'   lines(x, true, type = 'l', col = 2, lwd = 2, lty = 2)
#' }
#' legend("topleft", c("estimated", "truth"), col = 1:2, lty = 1:2, lwd = 2)
#' library(gamFactory)
#'
#' ###############
#' # [2] Example with single index
#' ###############
#' n <- 1000
#' 
#' X <- scale(matrix(rnorm(n * 5), n, 5))
#' alpha <- c(1, 0, 0, 0, 1)
#' si <- X %*% alpha
#' mu <- 10 * (0.5*cos(si) - 0.3*si)
#' 
#' x <- seq(-4, 4, length.out = n)
#' sigma =  .5+0.4*(x+4)*.5            # Scale
#' eps = 2*sin(x)                      # Skewness
#' del = 1 + 0.2*cos(3*x)              # Kurtosis
#' 
#' dat <-  mu + (del*sigma)*sinh((1/del)*asinh(qnorm(runif(n))) + (eps/del))
#' dataf <- data.frame(cbind(dat, x))
#' names(dataf) <- c("y", "x")
#' plot(x, dat, xlab = "x", ylab = "y")
#' 
#' dataf$X <- X
#' 
#' fit <- gam_nl(list(y ~ 
#'               s_nest(X, trans = trans_linear(), k = 20), ~ s(x), ~ s(x), ~ s(x)), 
#'               data = dataf, family = fam_shash())
#' 
#' plot(si, dataf$y)
#' lines(sort(si), mu[order(si)], col = 4)
#' points(si, as.matrix(fit$fitted.values)[ , 1], col = 2)
#' 
#' library(mgcViz)
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1)
#' plot(fit, inner = TRUE, select = 1)
#' 
#' fit <- getViz(fit, nsim = 50)
#' check(fit, type = "tnormal")
#' check0D(fit, type = "tnormal")
#' 
fam_shash <- function(){
  
  bundle <- bundle_shash()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}
