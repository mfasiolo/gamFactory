#'
#' Basic function for penalized regression
#' 
#' @description XXX.
#' @name penReg
#' @rdname penReg
#' @export penReg
penReg <- function(x, e, y) 
{
  if (sum(abs(e)) == 0) {
    b <- qr.coef(qr(x), y)
    b[!is.finite(b)] <- 0
    return(b)
  }
  qrx <- qr(x, LAPACK = TRUE)
  R <- qr.R(qrx)
  r <- ncol(R)
  rr <- Rrank(R)
  R[, qrx$pivot] <- R
  Qy <- qr.qty(qrx, y)[1:ncol(R)]
  k <- 0.01 * norm(R)/norm(e)
  qrr <- qr(rbind(R, e * k))
  edf <- sum(qr.Q(qrr)[1:r, ]^2)
  while (edf > 0.9 * rr) {
    k <- k * 10
    qrr <- qr(rbind(R, e * k))
    edf <- sum(qr.Q(qrr)[1:r, ]^2)
  }
  while (edf < 0.7 * rr) {
    k <- k/20
    qrr <- qr(rbind(R, e * k))
    edf <- sum(qr.Q(qrr)[1:r, ]^2)
  }
  b <- qr.coef(qrr, c(Qy, rep(0, nrow(e))))
  b[!is.finite(b)] <- 0
  b
}