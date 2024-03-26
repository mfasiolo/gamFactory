#'
#' Fit GAM models with non-linear effects
#' 
#' @param formula,family,data same arguments as in [mgcv::gam].
#' @param ... further arguments to be passed to [mgcv::gam].
#' @name gam_nl
#' @rdname gam_nl
#' @export
#'
gam_nl <- function(formula, family = fam_gaussian(), data = list(), fit = TRUE, ...){
  
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

