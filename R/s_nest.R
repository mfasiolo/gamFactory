#'
#' Defining nested smooths in GAM formulas
#' 
#' @name s_nest
#' 
#' @description Function used to define nested smooth effects. It works similarly to \link{mgcv::s}.
#'              
#' @param ... A \code{matrix} object containing variables that are used to define the nested smooth effect. 
#'            Each transformation has their specific requirements on how the variables in this matrix should be organized. 
#'            Use [data_linear], [data_mgks], [data_nexpsm] or [data_linear_nexpsm] to easily create the matrix.
#' @param trans The type of covariate transformation to use (see, [trans_linear], [trans_mgks], [trans_nexpsm], or [trans_linear_nexpsm]).
#' @param k The number of basis functions used to construct the outer smooth effect. 
#' @param m A vector of two positive integers. The first indicates the order of the 
#'          outer B-spline basis (see \link{basis_bspline}), the second is the order of the corresponding P-spline
#'          penalty. 
#' @param sarg a list of further arguments to be passed to the effect constructor. 
#' @rdname s_nest
#' @export
#' 
#' @details
#' The types of transformations currently provided are four: the linear transformation, the multivariate kernel 
#' smooth transformation, the exponential smoothing transformation and double nested effect 
#' which combines single index and exponential smooth. Specifically, 
#' \itemize{
#'  \item{\code{trans_linear}{ implements a linear transformation \eqn{X^\top \alpha}, the simplest type of
#'    transformation which can be used to specify a single index vector (i.e. a projection) or a linear effect.}}
#'  \item{\code{trans_mgks}{ implements a multivariate kernel smooth transformation based on the response variable
#'    observation vector \code{z} and corresponding distance matrix. The transformation is defined as
#'    \deqn{\tilde{s}(\mathbf{x}_i) = \frac{\sum_{j\in \mathcal{N}_i}K_{\mathbf{a}}(\mathbf{x}_{i},
#'    \mathbf{x}_{j})z_{j}}{\sum_{q\in \mathcal{N}_i}K_{\mathbf{a}}(\mathbf{x}_{i},\mathbf{x}_{q})},}
#'    with \eqn{K_{\mathbf{a}}} denoting the kernel of a multivariate p.d.f. parameterised by
#'    \eqn{\mathbf{a}}, and \eqn{\mathcal{N}_i} the index set of neighbours of \eqn{\mathbf{x}_i}.}}
#'  \item{\code{trans_nexpsm}{ implements an exponential smoothing transformation, defined as
#'    \deqn{\tilde{s}(x_i) = \tilde{s}_i = \omega_{i}\tilde{s}_{i-1}+(1-\omega_{i})x_{i}, \qquad \text{for } i \geq 1,}
#'    with \eqn{\tilde{s}_0 = x_0} and \eqn{\omega_i \in (0, 1)}. The smoothing factor can be defined as
#'    \eqn{\omega_i = \phi(\tilde{\mathbf{x}}_i^\top \mathbf{a})}, where \eqn{\phi} is the logistic function, 
#'    \eqn{\tilde{\mathbf{x}}_i} is a given vector.}}
#'  \item{\code{trans_linear_nexpsm}{ implements a combination of a linear transformation and an exponential
#'    smoothing transformation. The trasformation is defined as 
#'    \deqn{\tilde{s}({\bf x}_i, {\bf w}_i) = \tilde{s}_i = \omega_i \tilde{s}_{i-1} + (1-\omega_i)z_i,}
#'    where \eqn{z_i = {\bf x}_i^\top \bm a_{\mathrm{si}}}, for \eqn{i \geq 1,\ \tilde{s}_0 = z_0}, 
#'    and the smoothing factor is defined as \eqn{\omega_i = \phi({\bf w}_i^\top \bm a_{\mathrm{nexp}}),\ \omega_i \in (0,1)}. 
#'    As in the previous case, the smoothing factor can be defined as
#'    \eqn{\omega_i = \phi(\tilde{\mathbf{x}}_i^\top {\mathbf{a}_{\mathrm{si}})}}, where \eqn{\phi} is the logistic function
#'    and \eqn{\tilde{\mathbf{x}}_i} is a given vector.
#'    }}}
#' \code{trans_nexpsm} and \code{trans_linear_nexpsm} allow for different time granularities in the inner exponential and outer smooths. 
#' It is particularly useful for handling data structures where the response variable is observed at a 
#' different frequency than the covariates driving the smoothing rate (e.g., daily response with hourly covariates).
#' 
#' @examples
#' #####
#' # Example using trans_linear
#' ####
#' \dontrun{
#' library(gamFactory);library(mgcViz)
#' set.seed(6436)
#' n <- 1000; p <- 3; b <- 1:p
#' 
#' X <- matrix(runif(p * n), ncol = p)
#' y <- 2 * sin(X %*% b) + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1)
#' dat$X <- X # single-index predictors
#' 
#' # Fit the model
#' fit <- gam_nl(list(y~s_nest(X, trans = trans_linear(), k = 10, m = c(4,2)), ~1), 
#'               data = dat, family = fam_gaussian(), control=list(trace=TRUE))
#' 
#' # Plot the fit
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1) # plot smooth effects
#' print(plot(fit, inner = TRUE), pages = 1) # plot inner components
#' #####
#' # Example using trans_mgks
#' ####
#' n0 <- 50;n <- 1000
#' X0 <- cbind(runif(n0, -1, 1), runif(n0, -4, 4))
#' X <- cbind(runif(n, -1, 1), runif(n, -4, 4))
#' dist <- lapply(1:ncol(X), function(dd){
#'   t(sapply(1:nrow(X), function(ii) (X[ii, dd] - X0[ , dd])^2 ))
#' })
#' dist <- dist[[1]] + dist[[2]]
#' trueF <- function(x) 3 * x[ , 1] + x[ , 2]^2
#' temp <- t(sapply(1:n, function(ii) trueF(X0)))
#' 
#' z.mgks <- trueF(X)
#' dat <- data.frame(y = z.mgks*cos(z.mgks/2)+rnorm(n))
#' dat$Xks <- cbind(temp,dist)
#' colnames(dat$Xks) <- c(rep("y", nrow(X0)), rep("d1", nrow(X0)))
#' fit <- gam_nl(list(y~s_nest(Xks, trans = trans_mgks(), k = 10, m = c(4,2)), ~1), 
#'              data = dat, family = fam_gaussian(), control=list(trace=TRUE))
#' # Plot the fit
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1) # plot smooth effects
#' print(plot(fit, inner = TRUE), pages = 1) # plot inner components
#' #####
#' # Example using trans_nexpsm
#' ####
#' tim <- seq(0, 2*pi, length.out = n)
#' temp <- sin(tim) + rnorm(n, 0, 0.5)
#' Xi <- cbind(1, tim, tim^2)
#' Xi[,-1] <- scale(Xi[,-1], scale=FALSE)
#' aa.expsm <- c(2, 0.1, 0.1)
#' z.exp <- expsmooth(y = temp, Xi = Xi, beta = aa.expsm)$d0
#' dat <- data.frame(y = z.exp + sin(3 * z.exp)+rnorm(n))
#' dat$Xexp <- cbind(temp, Xi)
#' colnames(dat[,"Xexp"]) <- c("y", rep("x", ncol(Xi)))
#' 
#' fit <- gam_nl(list(y~s_nest(Xexp, trans = trans_nexpsm(), k = 10, m = c(4,2)), ~1), 
#'               data = dat, family = fam_gaussian(), control=list(trace=TRUE))
#' # Plot the fit
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1) # plot smooth effects
#' print(plot(fit, inner = TRUE), pages = 1) # plot inner components
#' #####
#' # Example using trans_linear_nexpsm
#' ####
#' n <- 1000; p <- 3; b_si <- 1:p; b_nexp <- 1:3
#' 
#' tim <- seq(0, 2*pi, length.out = n)
#' temp <- sin(tim) + rnorm(n, 0, 0.5)
#' w <- cbind(1, tim, tim^2)
#' w[,-1] <- scale(w[,-1], scale=FALSE)
#' x <- matrix(runif(p * n), ncol = p)
#' xw <- cbind(x,w)
#' colnames(xw) <- c("x1","x2","x3","w1","w2","w3")
#' y <- gamFactory:::deriv_si_nexp(X_si = x, X_nexp = w, param = c(b_si,b_nexp))$d0 + rnorm(n)
#' dat <- data.frame(y = y, xw = I(xw))
#' 
#' # Fit the model
#' fit <- gam_nl(list(y~s_nest(xw, trans = trans_linear_nexpsm(n_si = 3), k = 10, m = c(4,2)), ~1), 
#'               data = dat, family = fam_gaussian(), control=list(trace=TRUE))
#' 
#' # Plot the fit
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1) # plot smooth effects
#' print(plot(fit, inner = TRUE), pages = 1) # plot inner components
#' #####
#' }
#'
s_nest <- function(..., trans, k, m, sarg){
  
  vars <- as.list(substitute(list(...)))[-1]
  
  if( missing(sarg) ){ sarg <- list() }
  
  if( is.null(sarg$xt) ){ sarg$xt <- list() }
  sarg$xt$sumConv <- FALSE
  sarg$bs <- trans$type
  trans$type <- NULL
  if( length(trans) ){ sarg$xt$si <- trans }
  
  if( missing(k) ){
    k <- 10 + 7 # Because we impose 7 constraints
  } else {
    k <- k + 7
  }
  sarg$k <- k 
  
  if( missing(m) ){ m <- c(4, 2) }
  if( length(m) == 1 ){ m <- c(m, 2) }
  if( m[1] < 4 ){ message("m[1] should be >= 4, lower values might lead to convergence issues.") }
  sarg$m <- m
  
  # Call mgcv::s() to create smooth
  sm <- do.call("s", c(vars, sarg))
  
  # Create formula for term "s(...)" that can be used to build GAM model formula
  vnam <- as.character(vars)
  form <- do.call("call", c("s", c(vnam, sarg)))
  form <- paste(deparse(form), collapse="")
  
  # Transform "s(\"var1\", \"var2\", ...)" into "s(var1, var2, ...)" 
  for(ii in 1:length(vnam)){
    form <- gsub("\\s+", " ", gsub(paste0("\"",vnam[ii],"\""), vnam[ii], form, fixed=TRUE))
  }
  
  sm$form_term <- form
  
  return(sm)
  
}