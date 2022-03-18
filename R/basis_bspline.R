#'
#' Get B-spline design matrix constructor
#' 
#' @param k number of knots.
#' @param m order of the B-spline basis and penalty. For instance \code{m = 3} creates 
#'          a 3rd degree basis with a second order penalty.
#' @param lim range over which the basis is created.
#' 
#' @return A function that take arguments \code{x} and \code{deriv} and returns a P-Spline model matrix and its first
#'         \code{deriv} derivatives w.r.t. \code{x}.
#' 
#' @name basis_bspline
#' @rdname basis_bspline
#' @export basis_bspline
#' @examples
#' library(gamFactory)
#' 
#' tmp <- basis_bspline(k = 6, 
#'                      m = 3, 
#'                      lim = c(-1, 1)
#'                      )$eval(x = seq(-0.9, 0.9, 0.01), deriv = 3)
#'                      
#' par(mfrow = c(2, 2)) 
#' matplot(tmp$X0, type = 'l') # Design matrix
#' matplot(tmp$X1, type = 'l') # 1st derivative w.r.t. x 
#' matplot(tmp$X2, type = 'l') # 2nd derivative w.r.t. x 
#' matplot(tmp$X3, type = 'l') # 3rd derivative w.r.t. x
#' 
basis_bspline <- function(k, m, lim){
  
  force(k); force(m); force(lim);
  
  out <- list()
  
  .eval <- function(x, deriv){
    withCallingHandlers({
    gamFactory:::.basis_bspline(x = x, k = k, m = m, lim = lim, deriv = deriv)
    }, warning = function(w) {
      if (length(grep("there is \\*no\\* information about some basis coefficients", conditionMessage(w)))){
        invokeRestart("muffleWarning")
      }
    })
  }
  
  out <- structure(list("eval" = .eval), class = c("B-spline", "basis"))
  
  return(out)
}


##########################
# Internal to construct P-spline design
# 
.basis_bspline <- function(x, k, m, lim, deriv){
  
  n <- length(x)
  
  # Order of the penalty does not matter here (we just want design matrix)
  m <- c(m, 0)
  
  # Observations falling outside the boundary knots
  whIn <- x >= lim[1] & x <= lim[2]
  
  # Construct basis and penalty only using data inside knots
  sm <- smoothCon(object = s(x, bs = "ps", k = k, m = m), 
                  data = data.frame(x = x[whIn]), 
                  knots = list(x = lim), scale.penalty = FALSE)[[1]]
  
  # Get full design matrix using also data outside knots: need to call this to get X1, X2 and X3 
  X0 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T)$design
  
  if( any(abs(X0[whIn, ] - sm$X) > 1e-6)  ){ 
    stop("Problem in the creation of the P-spline design matrix") 
  }
  
  X1 <- X2 <- X3 <- NULL
  if(deriv > 0){
    X1 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T, derivs = x*0 + 1)$design
    if(deriv > 1){
      X2 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T, derivs = x*0 + 2)$design
    }
    if(deriv > 2){
      X3 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T, derivs = x*0 + 3)$design
    }
  }
  
  return( list("X0" = X0, "X1" = X1, "X2" = X2, "X3" = X3) )
}