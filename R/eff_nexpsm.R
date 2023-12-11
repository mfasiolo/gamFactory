#'
#' Build nested adaptive exponential smooth effect
#' 
#' @name eff_nexpsm
#' @rdname eff_nexpsm
#' @export eff_nexpsm
#'
eff_nexpsm <- function(y, Xi, basis, x0 = NULL, times = NULL){
  
  force(y); force(Xi); force(basis); force(x0); force(times);
  
  incall <- as.expression(quote(do.call("expsmooth", list("y" = y, "Xi" = Xi, "beta" = alpha, "times" = times, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_nexpsm", list("y" = y, "Xi" = Xi, "basis" = basis, "x0" = x0, "times" = times), quote = TRUE)))
  
  .eval <- .get_eff_eval_general()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("nexpsm", "nested"))
  
  return( out )
  
}









