#'
#' Basic function for penalized regression
#' 
#' @description XXX.
#' @name penReg
#' @rdname penReg
#' @export penReg
penReg <- function(x,e,y) {
  ## Taken from Simon N. Wood (in mgcv it is "pen.reg")
  ## get coefficients of penalized regression of y on matrix x
  ## where e is a square root penalty. Idea is to use e mainly for 
  ## regularization, so that edf is close to rank of x.  
  if (sum(abs(e))==0) { ## no penalization - easy
    b <- qr.coef(qr(x),y);b[!is.finite(b)] <- 0
    return(b)
  }
  ## need to adjust degree of penalization, so best to QR
  ## the x matrix up front...
  qrx <- qr(x,LAPACK=TRUE)
  R <- qr.R(qrx)
  r <- ncol(R)
  rr <- Rrank(R) ## rank of R/X
  R[,qrx$pivot] <- R ## unpivot
  Qy <- qr.qty(qrx,y)[1:ncol(R)]  
  ## now we want estimates with penalty weight low enough 
  ## EDF is k * rr where k is somewhere in e.g. (.7,.9)
  k <- .01 * norm(R)/norm(e)
  qrr <- qr(rbind(R,e*k));
  edf <- sum(qr.Q(qrr)[1:r,]^2) 

  rp <- Rrank(qr.R(qr(e))) ################## MATTEO Rank of penalty
  rov <- rp - (Rrank(qr.R(qrr)) - rr) ####### MATTEO number of dimensions in X that are penalized
  
  ################# MATTEO I added "&& edf > rr - 0.5*rov" so we stop in the "middle" of the range
  #################        of the penalty. Otherise, if 0.9*rr < rr - rov the while loop continue to 
  #################        run until we k = infinity !! The second while loop should be fine.
  while (edf > .9*rr && edf > rr - 0.5*rov) { ### increase penalization
    k <- k*10
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  while (edf<.7*rr) { ## reduce penalization
    k <- k/20
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  b <- qr.coef(qrr,c(Qy,rep(0,nrow(e))));b[!is.finite(b)] <- 0
  b
} ## pen.reg