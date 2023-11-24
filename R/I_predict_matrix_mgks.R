#'
#' Predict using MGKS smooth effects
#' 
#' @noRd
.predict.matrix.mgks <- function(object, data){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  if( is.null(si$xm) ){
    si$xm <- 0
  }
  
  X0 <- si$X0
  y0 <- si$y0
  Xi <- data[[object$term]]
  d <- ncol(X0)
  n <- nrow(Xi)
  n0 <- nrow(X0)
  
  # If TRUE then first n0 columns of Xi contain data y0 to be kernel smoother and
  # remaining "d" columns give the location at which we evaluate the kernel smooth.
  # If FALSE we only have the "d" columns and y0 was defined via "y0" argument in trans_mgks
  if( ncol(Xi) > d ){
    if( !is.null(y0) ){ stop(paste(object$term, "should have", d, "columns")) }
    y0 <- Xi[ , 1:(ncol(Xi)-d)]
    Xi <- Xi[ , -(1:(ncol(Xi)-d))]
    if(ncol(y0) != n0){
      stop(paste(object$term, "should have", n, "rows and", n0+d, "columns"))
    }
  } else {
    if( is.null("y0") ){ stop("Argument y0 missing in trans_mgks.") }
  }

  alpha <- si$alpha
  a0 <- alpha[1]
  a1 <- alpha[-1]
  xsm <- exp(a0) * mgks(y = y0, X = Xi, X0 = X0, beta = a1)$d0 
  
  # Compute outer model matrix
  X1 <- object$xt$basis$evalX(x = xsm - si$xm, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X1), length(alpha)), X1) 
  
  attr(Xtot, "inner_linpred") <- xsm
  
  return(Xtot)
  
}
