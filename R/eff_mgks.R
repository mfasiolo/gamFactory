#'
#' Build nested adaptive exponential smooth effect
#' 
#' @name eff_mgks
#' @rdname eff_mgks
#' @export eff_mgks
#'
eff_mgks <- function(y, dist, basis){
  
  force(y); force(dist); force(basis);
  
  na <- length(dist) + 1
  Xi <- NULL
  
  incall <- as.expression(quote(do.call("mgks", list("y" = y, "dist" = dist, "beta" = alpha, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_mgks", list("y" = y, "dist" = dist, "basis" = basis), quote = TRUE)))
  
  .eval <- .get_eff_eval_general()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("mgks", "nested"))
  
  return( out )
  
}

