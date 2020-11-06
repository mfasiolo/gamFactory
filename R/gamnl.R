#'
#' Fit GAM model with non-linear effects
#' 
#' @param formula,family,data same arguments as in [mgcv::gam].
#' @param ... further arguments to be passed to [mgcv::gam].
#' @name gamnl
#' @rdname gamnl
#' @export
#'
gamnl <- function(formula, family = gaulss(), data = list(), ...){
  
  out <- gam(formula = formula, family = family, data = data, fit = FALSE, ...)
  
  info <- .prepEffInfo(o = out)
  
  fam <- .buildGamlssSI(fam = bundle_gaussian(), effInfo = info)
  
  out$family <- fam()
  
  out <- gam(G = out, ...)

  out <- gamFactory:::.postProcNonLinGAM(o = out, info = info)
  
  return( out )
  
}

