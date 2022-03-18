#'
#' Build nested adaptive exponential smooth effect
#' 
#' @name eff_nexpsm
#' @rdname eff_nexpsm
#' @export eff_nexpsm
#'
eff_nexpsm <- function(y, Xi, basis, x0 = NULL){
  
  force(y); force(Xi); force(basis); force(x0)
  
  incall <- as.expression(quote(do.call("expsmooth", list("y" = y, "Xi" = Xi, "beta" = alpha, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_nexpsm", list("y" = y, "Xi" = Xi, "basis" = basis, "x0" = x0), quote = TRUE)))
  
  .eval <- .get_eff_eval_general()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("nexpsm", "nested"))
  
  return( out )
  
}









