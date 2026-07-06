#'
#' Log-likelihood of Generalized Pareto Distribution
#'
#' @description Log-likelihood of the Generalized Pareto Distribution (GPD), and its
#'              derivatives with respect to \code{param = cbind(phi, xi)}, where
#'              \code{xi} is the shape and \code{phi = sigma * (1 + xi)} is a
#'              reparametrisation of the scale \code{sigma} (chosen so that the
#'              parameters are close to orthogonal, see [`gamFactory::bundle_gpd`]).
#'              Returned in the list format used by [`gamFactory::llk_gaussian`]
#'              and friends.
#' @param y a vector of observations (non-negative).
#' @param param a matrix (or list) with 2 columns (elements), containing \code{phi} and
#'              \code{xi}, in this order.
#' @param deriv integer between 0 and 4 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-4
#'              additionally return \code{d1}-\code{d4}.
#' @name llk_gpd
#' @rdname llk_gpd
#' @export llk_gpd
#' @examples
#' library(gamFactory)
#' n <- 5
#' param <- c(1.5, 0.3) # phi and xi
#' y <- rexp(n)
#'
#' # Derivatives of GPD log-lik up to order 4
#' llk_gpd(y = y, param = param, deriv = 4)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_gpd(y = y, param = param, deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 1)$d1))
#'
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 2)$d2))
#'
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 3)$d3))
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = param, ord = 1:4)
#'
llk_gpd <- function(y, param, deriv = 0, ...) {
  
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")
  
  phi <- param[ , 1, drop = TRUE]
  xi <- param[ , 2, drop = TRUE]
  
  tol <- 1e-7
  if( any(abs(xi) < tol) ){ 
    xi[xi < tol & xi > 0] <- tol
    xi[xi > -tol & xi <= 0] <- - tol
  }
  
  Cxi <- xi * (1 + xi) * y
  
  out <- list()
  
  out[["d0"]] <- log1p(xi) - log(phi) - (1/xi+1) * log1p(Cxi / phi)
  
  if( deriv > 0 ){
    d1 <- (-phi+y+xi*y)/(phi*(phi+Cxi))
    d2 <- 1/(1+xi)+1/xi^2*log1p(Cxi/phi) - (1/xi+1)*(2*xi+1)*y/(phi+Cxi)
    out[["d1"]] <- list(d1, d2)
    
    if( deriv > 1 ){
      d11 <- (phi^2-2*phi*(1+xi)*y-xi*(1+xi)^2*y^2)/(phi^2*(phi+Cxi)^2)
      d12 <- -(1+xi)*y*(-2*phi+y+xi*y)/(phi*(phi+Cxi)^2)
      d22 <- - 1/(1+xi)^2-2/xi^3*log1p(Cxi/phi) + 
        2/xi^2*(2*xi+1)*y/(phi+Cxi) -
        (1/xi+1)*(2*y*(phi+Cxi)-(2*xi+1)^2*y^2)/(phi+Cxi)^2
      out[["d2"]] <- list(d11, d12, d22) 
      
      if( deriv > 2 ){
        A <- phi+xi*(1+xi)*y
        term_1 <- 2/(1+xi)^3-2*(1+xi)*(y+2*xi*y)^3/(xi*A^3)
        term_2 <- 6*(1+xi)*(1+2*xi)*y^2/(xi*A^2)-3*(y+2*xi*y)^2/(xi^2*A^2)
        term_3 <- 6*y/(xi^2*A)-6*(y+2*xi*y)/(xi^3*A)
        term_4 <- 6*log1p(Cxi/phi)/xi^4
        d111 <- (-2*phi^3+6*phi^2*(1+xi)*y+6*phi*Cxi^2/xi+2*Cxi^3/xi)/(phi^3*(phi+Cxi)^3)
        d112 <-(1+xi)*y*(-4*phi^2+3*phi*(1+xi)*y+xi*(1+xi)^2*y^2)/(phi^2*(phi+Cxi)^3)
        d122 <- 2*y*(phi^2-3*phi*(1+xi)^2*y+(1+xi)^3*y^2)/(phi*(phi+Cxi)^3)
        d222 <- term_1 + term_2 + term_3 + term_4
        out[["d3"]] <- list(d111, d112, d122, d222)

        if( deriv > 3 ){
          # AI GENERATED
          # 4th order derivatives, obtained via symbolic differentiation (sympy) of d0
          # w.r.t. phi and xi, then simplified via common-subexpression elimination.
          # Verified against finite differences of d3 (see llk_gpd's examples).
          g1 <- xi + 1
          g2 <- g1*xi
          g3 <- A^4
          g4 <- xi^3
          g5 <- g1^5
          g6 <- g5*y^4
          g7 <- g1^2
          g8 <- A^3
          g9 <- g8*y
          g10 <- g1^3
          g11 <- y^2
          g12 <- g11*xi
          g13 <- A^2
          g14 <- 6*g13
          g15 <- xi^2
          g16 <- y^3
          g17 <- g15*g16
          g18 <- g1^4
          g19 <- 4*A*g18
          g20 <- 1/g3
          g21 <- 6*g20
          g22 <- 2*xi
          g23 <- g22 + 1
          g24 <- 3*g10*g17
          g25 <- 3*xi
          g26 <- g13*y
          g27 <- g1*g26
          g28 <- A*g11*g22*g7
          g29 <- 2*g20
          g30 <- g29*y
          g31 <- 1/xi
          g32 <- g23^2
          g33 <- -g7
          g34 <- g1*g23
          g35 <- 2*g34
          g36 <- g2 - g22 + g35
          g37 <- 1/phi
          g38 <- g23^3
          g39 <- A*Cxi
          g40 <- xi^5
          g41 <- xi^4
          g42 <- 12*g18
          g43 <- g42*g9
          g44 <- g11*g14
          g45 <- g41*g5

          d1111 <- g21*(g10*g12*g14 - g17*g19 + g3 + g4*g6 - 4*g7*g9)/phi^4
          d1112 <- -g1*g30*(g23*g24 + 6*g27*(g25 + 1) - g28*(9*xi + 4) - 6*g8)/phi^3
          d1122 <- -g30*g31*(g22*g8 - g24*g32 - g26*xi*(g36 + 6*g7 - 2) + g27*(g33 - g35 + 5*xi + 3) +
                    g28*(g22*g23 + g32 + g36 - 1))/phi^2
          d1222 <- g11*g21*g31*g37*(-g12*g38*g7 - g13*g2 + g13*(g25 + g33 - g34 + 2) +
                    g23*g39*(4*xi + 3) + g32*g39)
          d2222 <- -g29*(12*A*g16*g32*g45 + g11*g13*g23*g4*g42 - g15*g18*g32*g44 + g15*g43 -
                    g16*g19*g38*g4 - 3*g23^4*g41*g6 - g23*g43*xi + 3*g3*g40 + g3*g42*log1p(Cxi/phi) -
                    g44*g45)/(g18*g40)

          out[["d4"]] <- list(d1111, d1112, d1122, d1222, d2222)
        }
      }
    }
  }
  
  return( out )
}


