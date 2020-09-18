#'
#' Costruct P-spline design with derivatives
#' 
#' @name PsplineDesign
#' @rdname PsplineDesign
#' @export PsplineDesign
#'
PsplineDesign <- function(x, k, m, lim, deriv){
  
  n <- length(x)
  
  # Observations falling outside the boundary knots
  whIn <- x >= lim[1] & x <= lim[2]
  
  # Construct basis and penalty only using data inside knots
  sm <- smoothCon(object = s(x, bs = "ps", k = k, m = m), 
                  data = data.frame(x = x[whIn]), 
                  knots = list(x = lim), scale.penalty = FALSE)[[1]]
  
  # Square root of penalty
  B <- .getBmatrix(P = sm$S[[1]], r = sm$rank)
  
  # Get full design matrix using also data outside knots 
  X0 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T)$design %*% B
  
  if( any(abs(X0[whIn, ] - sm$X %*% B) > 1e-6)  ){ 
    stop("Problem in the creation of the P-spline design matrix") 
  }
  
  X1 <- X2 <- X3 <- NULL
  if(deriv > 0){
    X1 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T, derivs = x*0 + 1)$design %*% B
    if(deriv > 1){
      X2 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T, derivs = x*0 + 2)$design %*% B
    }
    if(deriv > 2){
      X3 <- splines::spline.des(sm$knots, x = x, ord = sm$m[1] + 2, outer.ok = T, derivs = x*0 + 3)$design %*% B
    }
  }
  
  return( list("X0" = X0, "X1" = X1, "X2" = X2, "X3" = X3) )
}



