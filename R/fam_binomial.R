#'
#' The binomial family
#' 
#' @name fam_binomial
#' @rdname fam_binomial
#' @export fam_binomial
#' @examples 
#' library(gamFactory)
#' 
#' #### [1] Example with standard effects
#' dat <- gamSim(1,n=2000,dist="binary")
#' 
#' # use "cr" basis to save time, with 2000 data...
#' form <- y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
#'   s(x3,bs="cr")
#' 
#' # Create binomial family and fit it
#' fit <- gam_nl(list(form), data = dat, family = fam_binomial(n = 1))
#' plot(fit, pages = 1, scale = FALSE)
#' 
#' plot(dat$f, qlogis(fit$fitted.values))
#' abline(0, 1, col = 2)
#' 
#' fit_2 <- gam(list(form), data = dat, family = fam_binomial(n = 1),
#'              optimizer = c("outer", "bfgs"))
#' 
#' err <- max(abs(fit$fitted.values - drop(fit_2$fitted.values)))
#' if(err > 1e-6){
#'   stop("Discrepancy between gam and gam_nl")
#' }
#' 
#' #### [2] Example with single index
#' n <- 2000
#' m <- 50
#' X <- scale(matrix(rnorm(n * 5), n, 5))
#' alpha <- c(1, 0, 0, 0, 1)
#' si <- X %*% alpha
#' eta <- 0.5*cos(si) - 0.3*si
#' y <- rbinom(n = n, size = m, p = plogis(eta))
#' plot(si, y)
#' 
#' my_data <- data.frame(y = y, si = si)
#' my_data$X <- X
#' 
#' fit <- gam_nl(list(y ~ s_nest(X, trans = trans_linear(), k = 20)), data = my_data, family = fam_binomial(n = m))
#' 
#' plot(si, y)
#' lines(sort(si), m*plogis(eta)[order(si)], col = 4)
#' points(si, m*as.matrix(fit$fitted.values)[ , 1], col = 2)
#' 
#' library(mgcViz)
#' fit <- getViz(fit)
#' plot(fit) + l_fitLine()
#' plot(fit, inner = TRUE)
#' 
fam_binomial <- function(n){
  
  bundle <- bundle_binomial(n)
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}

