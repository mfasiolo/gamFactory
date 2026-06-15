#' Build Nested Adaptive Exponential Smoothing Effect Evaluator
#' 
#' @description This function acts as a constructor for dual-nested smooth effects. 
#' It bundles the covariate data for both single-index projection and dynamic smoothing 
#' weights, along with initialization states and a basis expansion method, into an 
#' object that can evaluate the smoothed effect and its analytical derivatives dynamically.
#' 
#' @param X_si A numeric design matrix of dimensions \code{n} by \code{n_si} containing 
#' the covariates to be projected down into the single-index dimension.
#' @param X_nexp A numeric design matrix of dimensions \code{n} by \code{n_nexp} containing 
#' the covariates used to dynamically generate the adaptive exponential smoothing weights.
#' @param basis A structure or list containing a basis evaluation function 
#' \code{evalX(x, deriv)}. This function must accept a vector of smoothed values 
#' and return a list containing basis matrices (\code{X0}, \code{X1}, \code{X2}, \code{X3}) 
#' corresponding to the specified derivative order.
#' @param Z0 A numeric scalar or vector providing the initial state value used to 
#' start the exponential smoothing sequence.
#' @param times An optional numeric vector indicating specific time steps or indices 
#' for the sequential evaluation. Defaults to \code{NULL}.
#' @param alpha_center An optional numeric offset vector of length \code{n_si} added directly 
#' to the single index coefficients \code{alpha_si}. If \code{NULL}, no offset is applied.
#' @param positive_si A logical flag. If \code{TRUE}, it strictly enforces positivity 
#' on the single-index weights during dynamic evaluation. Defaults to \code{FALSE}.
#' 
#' @details The function maps covariates through a sequential three-layer transformation. 
#' The matrices \code{X_si} and \code{X_nexp} are combined with their respective parameters 
#' to construct a dynamically smoothed sequence, which is subsequently transformed via 
#' an outer smooth basis expansion.
#' 
#' For a given combined parameter vector \code{param = c(alpha_nexp, alpha_si, beta)}:
#' \enumerate{
#'   \item The single-index projection step is calculated as:
#'         \deqn{z_t = \mathbf{X}_{si, t \cdot} (\bm{\alpha}_{si} + \bm{\alpha}_{center})}
#'   \item The dynamic weights and the adaptive exponential smoothing sequence are calculated recursively as:
#'         \deqn{\omega_t = \text{sigmoid}(\mathbf{X}_{nexp, t \cdot} \bm{\alpha}_{nexp})}
#'         \deqn{s_t = \omega_t s_{t-1} + (1 - \omega_t) z_t, \qquad \text{with } s_0 = Z_0}
#'   \item The outer smooth function value is calculated via linear combination with the basis matrix evaluated at the smoothed sequence:
#'         \deqn{f(s_t) = \mathbf{B}_0(s_t) \bm{\beta}}
#' }
#' 
#' When \code{deriv >= 1}, partial derivatives of the spline basis with respect to 
#' both sets of inner projection parameters (\eqn{\bm{\alpha}_{nexp}} and \eqn{\bm{\alpha}_{si}}) 
#' are computed analytically via the chain rule and recursive differentiation across 
#' the time steps, yielding corresponding \eqn{f^{(m)}} components.
#' 
#' @return An object of class \code{c("si_nexpsm", "nested")} containing:
#' \itemize{
#'   \item{\code{eval(param, deriv = 0)}}{ A function designed to accept a combined 
#'         numeric parameter vector \code{param}. The first \code{ncol(X_nexp)} elements 
#'         represent the smoothing weight parameters \code{alpha_nexp}, the subsequent 
#'         \code{ncol(X_si)} elements represent the single-index weights \code{alpha_si}, 
#'         and the remainder represent the spline coefficients \code{beta}. It returns 
#'         a fully evaluated state object.}
#' }
#' @seealso [trans_linear_nexpsm], [s_nest]
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