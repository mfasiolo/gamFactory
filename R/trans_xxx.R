#'
#' Specifying inner transformations
#' 
#' @description Functions to specify nested covariate transformations, to be used
#'              within \link{s_nest}.

#' @rdname trans_xxx
#' 
#' @param pord The order of the differences of the transformation parameters
#'             that will be penalised. E.g., \code{pord = 1} corresponds to penalising the squared 
#'             first order differences \eqn{(alpha_k - alpha_{k-1})^2}. Ignored if argument \code{S}
#'             is supplied. 
#' @param S Positive semi-definite matrix used to specify a generalised ridge
#'          penalty on the transformation parameters. In \code{trans_linear}, if \code{S} is not supplied, the \code{pord}-order differences 
#'          between parameters are penalised, see the \code{pord} argument.
#' @param alpha Vector containing the initial values for the parameters of the transformation. If provided they can be used
#'              by the optimiser.
#' @param y0 Vector of \eqn{n} observations corresponding to the rows of the variables in \code{s_nest}.
#' 
#' @details The types of transformations currently provided are:
#' \itemize{
#'  \item{\code{trans_linear}{ a linear transformation \eqn{X^\top \alpha}, which can be used to specify, e.g, a single index vector (i.e. a projection) or
#'                       a linear effect.}}
#'  \item{\code{trans_mgks}{ a multivariate kernel smooth transformation based on the response variable observation 
#'                          vector \code{y0} and corresponding distance matrix.}}
#'  \item{\code{trans_nexpsm}{ an exponential smoothing transformation.}}
#' }
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

