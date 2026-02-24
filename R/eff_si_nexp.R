#'
#' Build nested adaptive linear transform and exponential smooth effect
#' 
#' @name eff_si_nexp
#' @rdname eff_si_nexp
#' @export eff_si_nexp
#'

eff_si_nexp <- function(X_si, X_nexp, basis, Z0, times = NULL, alpha_center = NULL, positive_si = FALSE) {
  force(X_si); force(X_nexp); force(basis); force(Z0); force(times); force(alpha_center); force(positive_si)
  # z0 is the initial value used to smooth the sequence. 
  # It is currently set to NULL by default, waiting for optimization later.
  
  n_si <- ncol(X_si)
  n_nexp <- ncol(X_nexp)
  na <- n_si + n_nexp  # length of alpha without alpha_scale

  incall <- as.expression(quote(
    do.call("deriv_si_nexp",
            list(
              "X_si"         = X_si,
              "X_nexp"       = X_nexp,
              "param"        = c(alpha_nexp, alpha_si),
              "times"        = times,
              "deriv"        = deriv,
              "alpha_center" = alpha_center,
              "Z0"           = Z0,          
              "positive_si"  = positive_si
            ),
            quote = TRUE)
  ))

  efcall <- as.expression(quote(
    do.call("eff_si_nexp",
            list(
              "X_si"         = X_si,
              "X_nexp"       = X_nexp,
              "basis"        = basis,
              "Z0"           = Z0,
              "times"        = times,
              "alpha_center" = alpha_center,
              "positive_si"  = positive_si
            ),
            quote = TRUE)
  ))
  
  
  .eval <- .get_eff_eval_si_nexp()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("si_nexpsm", "nested"))
  
  return( out )
  
}