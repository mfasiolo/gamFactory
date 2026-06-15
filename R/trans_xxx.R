#' Specifying inner transformations
#' 
#' @description Functions to specify nested covariate transformations, to be used
#'              within \link{s_nest}.
#'
#' @rdname trans_xxx
#' 
#' @param pord The order of the differences of the linear transformation parameters
#'             that will be penalised. E.g., \code{pord = 1} corresponds to penalising the squared 
#'             first order differences \eqn{(\alpha_k - \alpha_{k-1})^2}. Ignored if argument \code{S}
#'             is supplied. 
#' @param S Positive semi-definite matrix used to specify a generalised ridge
#'          penalty on the transformation parameters. In \code{trans_linear}, if \code{S} is not supplied, the \code{pord}-order differences 
#'          between parameters are penalised, see the \code{pord} argument.
#' @param S_si Penalty matrix used in trans_linear_nexpsm. For the linear transformation parameters.
#' @param S_nexp Penalty matrix used in expsmooth. For the linear transformation parameters.
#' @param alpha Vector containing the initial values for the parameters of the transformation.
#' @param alpha_si Initial values used in trans_linear_nexpsm. For linear transformation parameters.
#' @param alpha_nexp Initial values used in trans_linear_nexpsm. For expsmooth parameters.
#' @param y0 Vector of \eqn{n} observations corresponding to the rows of the variables in \code{s_nest}.
#' @param n_si Length of \code{alpha_si} if \code{alpha_si} is not provided.
#' @param n_nexp Length of \code{alpha_nexp} if \code{alpha_nexp} is not provided.
#' 
#' @details The types of transformations currently provided are four: the linear transformation, 
#' the multivariate kernel smooth transformation, the exponential smoothing transformation, and a double 
#' nested effect which combines a single index and exponential smooth. Specifically:
#' \itemize{
#'   \item{\code{trans_linear}}{ implements a linear transformation \eqn{X^\top \alpha}, the simplest type of
#'         transformation which can be used to specify a single index vector (i.e. a projection) or a linear effect.}
#'   \item{\code{trans_mgks}}{ implements a multivariate kernel smooth transformation based on the response variable
#'         observation vector \code{y0} and corresponding distance matrix. The transformation is defined as
#'         \deqn{\tilde{s}(\mathbf{x}_i) = \frac{\sum_{j\in \mathcal{N}_i}K_{\mathbf{a}}(\mathbf{x}_{i},
#'         \mathbf{x}_{j})y0_{j}}{\sum_{q\in \mathcal{N}_i}K_{\mathbf{a}}(\mathbf{x}_{i},\mathbf{x}_{q})},}
#'         with \eqn{K_{\mathbf{a}}} denoting the kernel of a multivariate p.d.f. parameterised by
#'         \eqn{\mathbf{a}}, and \eqn{\mathcal{N}_i} the index set of neighbours of \eqn{\mathbf{x}_i}.}
#'   \item{\code{trans_nexpsm}}{ implements an exponential smoothing transformation, defined as
#'         \deqn{\tilde{s}(x_i) = \tilde{s}_i = \omega_{i}\tilde{s}_{i-1}+(1-\omega_{i})x_{i}, \qquad \text{for } i \geq 1,}
#'         with \eqn{\tilde{s}_0 = x_0} and \eqn{\omega_i \in (0, 1)}. The smoothing factor can be defined as
#'         \eqn{\omega_i = \phi(\tilde{\mathbf{x}}_i^\top \mathbf{a})}, where \eqn{\phi} is the logistic function, 
#'         and \eqn{\tilde{\mathbf{x}}_i} is a given vector.}
#'   \item{\code{trans_linear_nexpsm}}{ implements a combination of a linear transformation and an exponential
#'         smoothing transformation. The transformation is defined as 
#'         \deqn{\tilde{s}({\bf x}_i, {\bf w}_i) = \tilde{s}_i = \omega_i \tilde{s}_{i-1} + (1-\omega_i)z_i,}
#'         where \eqn{z_i = {\bf x}_i^\top \bm a_{\mathrm{si}}}, for \eqn{i \geq 1,\ \tilde{s}_0 = z_0}, 
#'         and the smoothing factor is defined as \eqn{\omega_i = \phi({\bf w}_i^\top \bm a_{\mathrm{nexp}}),\ \omega_i \in (0,1)}. 
#'         As in the previous case, the smoothing factor can be defined as
#'         \eqn{\omega_i = \phi(\tilde{\mathbf{x}}_i^\top {\mathbf{a}_{\mathrm{si}})}}, where \eqn{\phi} is the logistic function
#'         and \eqn{\tilde{\mathbf{x}}_i} is a given vector.}
#' }
#' 
#' \code{trans_nexpsm} and \code{trans_linear_nexpsm} allow for different time granularities in the inner exponential and outer smooths. 
#' It is particularly useful for handling data structures where the response variable is observed at a 
#' different frequency than the covariates driving the smoothing rate (e.g., daily response with hourly covariates).
#' 
#' @export trans_linear
#'
trans_linear <- function(pord, S, alpha, a0){
  
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "si"
  
  return(out)
  
}

#
# Specifying a MGKS transformation
# 
#' @rdname trans_xxx
#' @export trans_mgks
#'
trans_mgks <- function(y0, alpha){
  
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "mgks"
  
  return(out)
  
} 
# trans_mgks <- function(X0, y0, alpha){
#  
#   if( missing(X0) ){ stop("Argument \"X0\" is missing but a value is required") }
#    
#   out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
#   out$type <- "mgks"
#  
#   return(out)
#   
# } 

#
# Specifying exponential smooth transformation 
# 
#' @rdname trans_xxx
#' @export trans_nexpsm
#'
trans_nexpsm <- function(S, alpha){
  
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "nexpsm"
  
  return(out)
  
} 

#
# Specifying linear transform + exponential smooth transformation 
# 
#' @rdname trans_xxx
#' @export trans_linear_nexpsm
#'
trans_linear_nexpsm <- function(
    pord         = NULL,
    S_si         = NULL,
    S_nexp       = NULL,
    alpha_nexp   = NULL,
    alpha_si     = NULL,
    alpha_scale  = NULL,
    center       = FALSE,
    alpha_center = NULL,
    n_si         = NULL,
    n_nexp       = NULL,
    Z0           = NULL,   # initialization for expsmooth
    positive_si = FALSE
){
  # if missing n_si/n_nexp，get it from alpha_si/alpha_nexp
  if (is.null(n_si)   && !is.null(alpha_si))   n_si   <- length(alpha_si)
  if (is.null(n_nexp) && !is.null(alpha_nexp)) n_nexp <- length(alpha_nexp)
  if (is.null(alpha_nexp) && is.null(alpha_si) && is.null(n_si) && is.null(n_nexp)) {
    stop("You must provide at least one of: alpha_si, alpha_nexp, n_si, n_nexp.")
  }
  
  out <- as.list(environment())
  out$type <- "si_nexpsm"
  return(out)
}








