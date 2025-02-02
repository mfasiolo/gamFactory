#'
#' Fit GAM models with non-linear effects
#' 
#' @description This function fits a generalized additive model for location scale and shape (GAMLSS) 
#'              with non-linear effects using [mgcv::gam]
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
#'               data = dat, family = fam_gaussian(), control=list(trace=TRUE))
#' 
#' # Plot the fit
#' fit <- getViz(fit)
#' print(plot(fit), pages = 1) # plot smooth effects
#' print(plot(fit, inner = TRUE), pages = 1) # plot inner components
#'
#'
gam_nl <- function(formula, family = fam_gaussian(), data = list(), fit = TRUE, ...){
  
  if( !is.list(formula) ){
    formula <- list(formula)
  }
  
  ddd <- match.call(expand.dots = FALSE)$`...`
  
  out <- ddd$G
  if( is.null(out) ){ # Do not build if G already provided
    form_comp <- .compile_formula(formula)
    
    out <- gam(formula = form_comp, family = family, data = data, fit = FALSE, ...)
    
    info <- prep_info(o = out)
    
    fam <- build_family_nl(bundle = do.call(family$bundle_nam, as.list(family$store)), info = info)
    
    out$family <- fam()
  }
  
  if( fit ){
    out <- gam(G = out, ...)

    out <- postproc_gam_nl(o = out, info = info)
    
    for(ii in 1:length(formula)){ # because form_comp evaluates all the arguments
      out$formula[[ii]] <- formula[[ii]]
    }
  }
  
  class(out) <- c("gamnl", class(out))
  
  return( out )
  
}

