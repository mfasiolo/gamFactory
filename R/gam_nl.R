#'
#' Fit GAM model with non-linear effects
#' 
#' @param formula,family,data same arguments as in [mgcv::gam].
#' @param ... further arguments to be passed to [mgcv::gam].
#' @name gam_nl
#' @rdname gam_nl
#' @export
#'
gam_nl <- function(formula, family = gaulss(), data = list(), ...){
  
  out <- gam(formula = formula, family = family, data = data, fit = FALSE, ...)
  
  info <- prep_info(o = out)
  
  fam <- build_family_nl(bundle = bundle_gaussian(), info = info)
  
  out$family <- fam()
  
  out <- gam(G = out, ...)

  out <- postproc_gam_nl(o = out, info = info)
  
  return( out )
  
}

