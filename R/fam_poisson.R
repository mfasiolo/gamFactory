#'
#' The Poisson family
#' 
#' @name fam_poisson
#' @rdname fam_poisson
#' @export fam_poisson
#' @examples 
#' library(gamFactory)
#' 
#' #### [1] Example with standard effects
#' dat <- gamSim(1,n=2000,dist="poisson",scale=.1)
#' 
#' # use "cr" basis to save time, with 2000 data...
#' form <- y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
#'   s(x3,bs="cr")
#' 
#' # Create poisson family and fit it
#' fit <- gam_nl(list(form), data = dat, family = fam_poisson())
#' plot(fit, pages = 1, scale = FALSE)
#' 
#' fit_2 <- gam(list(form), data = dat, family = fam_poisson(),
#'              optimizer = c("outer", "bfgs"))
#' 
#' err <- max(abs(fit$fitted.values - drop(fit_2$fitted.values)))
#' if(err > 1e-6){
#'   stop("Discrepancy between gam and gam_nl")
#' }
#' 
#' 
#' #### [2] Example with single index
#' n <- 2000
#' X <- scale(matrix(rnorm(n * 5), n, 5))
#' alpha <- c(1, 0, 0, 0, 1)
#' si <- X %*% alpha
#' lrate <- log(10*exp(0.5*cos(si) - 0.2*si))
#' y <- rpois(n = n, lambda = exp(lrate))
#' plot(si, y)
#' 
#' my_data <- data.frame(y = y, si = si)
#' my_data$X <- X
#' 
#' fit <- gam_nl(list(y ~ s_nest(X, trans = trans_linear(), k = 20)), data = my_data, family = fam_poisson())
#' 
#' plot(si, y)
#' lines(sort(si), exp(lrate)[order(si)], col = 4)
#' points(si, as.matrix(fit$fitted.values)[ , 1], col = 2)
#' 
#' library(mgcViz)
#' fit <- getViz(fit)
#' plot(fit) + l_fitLine()
#' plot(fit, inner = TRUE)
#' 
#' fit <- getViz(fit, nsim = 50)
#' check(fit)
#' check0D(fit, type = "deviance")
#' 
fam_poisson <- function(){
  
  bundle <- bundle_poisson()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}

