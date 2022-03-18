#'
#' Build nested adaptive exponential smooth effect
#' 
#' @name eff_mgks
#' @rdname eff_mgks
#' @export eff_mgks
#'
eff_mgks <- function(y, X, Xi, basis){
  
  force(y); force(X); force(Xi); force(basis);
  
  incall <- as.expression(quote(do.call("mgks", list("y" = y, "X" = X, "X0" = Xi, "beta" = alpha, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_mgks", list("y" = y, "X" = X, "Xi" = Xi, "basis" = basis), quote = TRUE)))
  
  .eval <- .get_eff_eval_general()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("mgks", "nested"))
  
  return( out )
  
}

