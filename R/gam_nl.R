#'
#' Fit GAM models with nested effects
#' 
#' @description This function fits a generalized additive model for location scale and shape (GAMLSS) 
#'              with nested effects using [mgcv::gam]
#' 
#' 
#' @param formula A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param family A family object specifying the distribution and link to use in fitting. 
#'                Available families are Gaussian ([fam_gaussian]), binomial ([fam_binomial]), 
#'                generalised Pareto distribution ([fam_gpd]), Poisson ([fam_poisson]) and 
#'                Sinh-Arsinh ([fam_shash]).
#' @param data A data frame or list that includes the model's response variable along with the covariates specified in the formula.
#'             For the structure required by the [s_nest] effects, refer to [trans_linear], [trans_mgks], or [trans_nexpsm]
#' @param fit Same argument as in [mgcv::gam]. If this argument is \code{TRUE} then \code{gam_nl} sets up the
#'            model and fits it, but if it is \code{FALSE} then the model is set up and an object containing what would be
#'            required to fit is returned.
#' @param n_init number of candidate directions tried when initialising each true
#'               single-index (\code{trans_linear}) nested effect. See
#'               \code{\link{build_family_nl}}.
#' @param n_eigen number of eigenvector-based candidate directions (out of \code{n_init})
#'                used when initialising each single-index nested effect. See
#'                \code{\link{build_family_nl}}.
#' @param oversample oversampling factor used to spread the remaining candidate
#'                   directions over the sphere. See \code{\link{build_family_nl}}.
#' @param ... further arguments to be passed to [mgcv::gam].
#' @name gam_nl
#' @rdname gam_nl
#' @return A \code{gamnl} and \code{gamObject}. See \code{?gamObject}.
#' 
#' @export
#'
#' @examples
#
#' #####
#' # Single index example
#' ####
#' library(gamFactory);library(mgcViz)
#' set.seed(6436)
#' n <- 1000; p <- 3; b <- 1:p
#' 
#' X <- matrix(runif(p * n), ncol = p)
#' x1 <- rnorm(n)
#' y <- 2 * sin(X %*% b) + x1^2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1)
#' dat$X <- X # single-index predictors
#' 
#' # Fit the model
#' fit <- gam_nl(list(y~s_nest(X, trans = trans_linear()) + s(x1), ~1),
#'               data = dat, family = fam_gaussian())
#' 
#' # Plot the smooth and single index vector
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1) # plot smooth effects
#' print(plot(fit, inner = TRUE), pages = 1) # plot inner components
#' 
#' # Plot the fit to data vs the true projected data
#' Xb <- dat$X %*% b
#' tmp_dat <- dat
#' tmp_dat$x1 <- 0
#' plot(Xb, dat$y)
#' lines(sort(Xb), predict(fit, newdata = tmp_dat)[order(Xb), 1], col = 2, lwd = 2)
#' lines(sort(Xb), 2 * sin(sort(Xb)), col = 4, lwd = 2)
#' 
#' #####
#' # Single index example (II)
#' # Same as before but X elements now follow a normal distribution, hence 
#' # projected data is more spread out and we need more basis functions.
#' 
#' X <- matrix(rnorm(p * n), ncol = p)
#' x1 <- rnorm(n)
#' y <- 2 * sin(X %*% b) + x1^2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1)
#' dat$X <- X # single-index predictors
#' 
#' # Fit the model, note k = 20
#' fit <- gam_nl(list(y ~ s_nest(X, k = 20, trans = trans_linear()) + s(x1), ~1),
#'               data = dat, family = fam_gaussian())
#' 
#' # Plot the fit
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1)
#' print(plot(fit, inner = TRUE), pages = 1)
#' 
#' # Plot the fit to data vs the true projected data
#' Xb <- dat$X %*% b
#' tmp_dat <- dat
#' tmp_dat$x1 <- 0
#' plot(Xb, dat$y)
#' lines(sort(Xb), predict(fit, newdata = tmp_dat)[order(Xb), 1], col = 2, lwd = 2)
#' lines(sort(Xb), 2 * sin(sort(Xb)), col = 4, lwd = 2)
#'
gam_nl <- function(formula, family = fam_gaussian(), data = list(), fit = TRUE,
                    n_init = 1000, n_eigen = 10, oversample = 10, ...){
  
  if( !is.list(formula) ){
    formula <- list(formula)
  }
  
  dots <- list(...)
  
  out <- dots$G
  
  ## Arguments for the first gam() call:
  ## remove sp so it is not used when fit = FALSE
  build_dots <- dots
  build_dots$sp <- NULL
  
  ## Arguments for the second gam() call:
  ## keep sp, but remove G because we pass G explicitly
  fit_dots <- dots
  fit_dots$G <- NULL
  fit_dots$data <- NULL
  
  if( is.null(out) ){ # Do not build if G already provided
    form_comp <- .compile_formula(formula)
    
    out <- do.call("gam", c(list(formula = form_comp, family = quote(family), data = quote(data), fit = FALSE), build_dots))
    
    info <- prep_info(o = out)
    
    fam <- build_family_nl(bundle = do.call(family$bundle_nam, as.list(family$store)), info = info, link = family$link,
                            n_init = n_init, n_eigen = n_eigen, oversample = oversample)
    
    out$family <- fam()
  }
  
  if( fit ){

    if( is.null(fit_dots$in.out) && is.null(dots$sp) ){
      fit_dots$in.out <- list("sp" = .nested_initial_sp(out), "scale" = 1)
    }

    out <- do.call("gam", c(list(G = out), fit_dots))

    out <- postproc_gam_nl(o = out, info = info)
    
    for(ii in 1:length(formula)){ # because form_comp evaluates all the arguments
      out$formula[[ii]] <- formula[[ii]]
    }
  }
  
  out$call <- match.call()
  
  class(out) <- c("gamnl", class(out))
  
  return( out )
  
}

