#'
#' Internal function for penalized regression
#' 
#' @description This is an internal function, which need to be exported for reasons having to 
#'              do with R environments. Please do not use it.
#' @param x design matrix
#' @param e matrix square root of penalty
#' @param y response vector
#' @name penreg
#' @rdname penreg
#' @export penreg
penreg <- function(x,e,y) {
  ## Taken from Simon N. Wood (in mgcv it is "pen.reg")
  ## get coefficients of penalized regression of y on matrix x
  ## where e is a square root penalty. Idea is to use e mainly for 
  ## regularization, so that edf is close to rank of x.  
  if (sum(abs(e))==0) { ## no penalization - easy
    b <- qr.coef(qr(x),y)
    b[!is.finite(b)] <- 0
    return(b)
  }
  ## need to adjust degree of penalization, so best to QR
  ## the x matrix up front...
  qrx <- qr(x, LAPACK=TRUE)
  R <- qr.R(qrx)
  r <- ncol(R)
  rr <- Rrank(R) ## rank of R/X
  R[ , qrx$pivot] <- R ## unpivot
  Qy <- qr.qty(qrx, y)[1:ncol(R)]  
  ## now we want estimates with penalty weight low enough 
  ## EDF is k * rr where k is somewhere in e.g. (.7,.9)
  k <- .01 * norm(R)/norm(e)
  qrr <- qr(rbind(R, e*k));
  edf <- sum(qr.Q(qrr)[1:r, ]^2)
  ## compute rank of e less rank of the space penalized by e not in
  ## range space of x, this is how much penalty can in principle change
  ## edf by. Needed for corner cases where e.g. penalty is imposed for
  ## identifiability reasons and then only penalizes null space of x...
  re <- min(sum(colSums(abs(e))!=0),nrow(e)) - Rrank(qr.R(qrr)) + rr
  while (edf > rr-.1*re) { ## increase penalization
    k <- k*10
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  while (edf<0.85*rr) { ## reduce penalization (was 0.7)
    k <- k/5 ## was 20! 
    qrr <- qr(rbind(R, e*k));
    edf <- sum(qr.Q(qrr)[1:r, ]^2)
  } 
  b <- qr.coef(qrr,c(Qy,rep(0,nrow(e))))
  b[!is.finite(b)] <- 0
  return(b)
} ## pen.reg