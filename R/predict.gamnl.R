#'
#' Predictions from GAM models with non-linear effects
#' 
#' @param formula,family,data same arguments as in [mgcv::gam].
#' @param ... further arguments to be passed to [mgcv::gam].
#' @name predict.gamnl
#' @rdname predict.gamnl
#' @export
#'
predict.gamnl <- function(object, newdata, type="link", se.fit=FALSE, ...){
  
  out <- predict.gam(object = object, newdata = newdata, type = type, block.size = Inf, ...)
  
  return(out)
  
}