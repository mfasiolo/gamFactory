#'
#' Defining nested smooths in GAM formulae
#' 
#' @name s_nest
#' 
#' @description Function used to define nested smooth effects. It works similarly to \link{mgcv::s}.
#'              
#' @param ... A \code{matrix} object containing variables that are used to define the nested smooth effect. See Details
#' @param trans The type of covariate transformation to use (see, [trans_linear], [trans_mgks], or [trans_nexpsm]).
#' @param k The number of basis functions used to construct the outer smooth effect. 
#' @param m A vector of two positive integers. The first indicates the order of the 
#'          outer B-spline basis (see \link{basis_bspline}), the second is the order of the corresponding P-spline
#'          penalty. 
#' @param sarg a list of further arguments to be passed to the effect constructor. 
#' @rdname s_nest
#' @export
#' 
#' @details
#' The types of transformations currently provided are three: the linear transformation, the multivariate kernel 
#' smooth transformation and the exponential smoothing transformation. Here we provide details on the definition 
#' of the list of variables used to define each nested effect. In the following we will refer to the list 
#' of variables passed as \code{X}.
#' 
#' If \code{trans = trans_linear()}, \code{X} is a \eqn{n\times p} matrix with a covariate for each column.
#' 
#' When \code{trans = trans_mgks()}, the matrix \code{X} is organized into two distinct blocks, \code{X = cbind(Y, D)}.  
#' The \eqn{j}-th column in \code{Y} corresponds directly to the \eqn{j}-th column in \code{D}, with the two 
#' blocks defined as follows:  
#'  
#' - The columns in \code{Y} are all named with \code{"y"} and represent observations of the variable that is first extrapolated and then smoothed.  
#' - The columns in \code{D} are named with \code{"d"} prefixes, and the \eqn{j}-th column in \code{D} contains the distance between the location where the
#'  corresponding \eqn{j}-th column in \code{Y} was observed and the location where all other covariates in the model were observed.  
#' 
#' When \code{trans = trans_nexpsm()}, the matrix \code{X} is divided into three distinct blocks: \code{X = cbind(y, Xi, times)}. The third block, \code{times}, is optional.  
#'  
#' - The columns in \code{y} are all named with \code{"y"} and represent the variable to be exponentially smoothed (stored as a vector or matrix). These columns may have a higher resolution compared to other covariates (e.g., hourly resolution for \code{y} versus daily resolution for other covariates and the response variable).  
#' - The columns in \code{Xi} are all named with \code{"x"} and correspond to the model matrix used to estimate the exponential smoothing rate.  
#' - The column in \code{times}, named \code{"times"}, contains the vector of times at which the response variable of the GAM is observed.  
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
#' 
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