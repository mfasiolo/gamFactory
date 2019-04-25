#'
#' Calculating derivatives of log-likelihood wrt regression coefficients
#'
#' @description Intended for use with stack effects. 
#' Converts derivatives of the log-likelihood with respect to
#' the non linear predictor eta into derivatives with respect to the
#' regression coefficients and smoothing parameters
#'
#' @param beta regression coefficients
#' @param theta vector of other parameters of the response variable distribution (to be included in beta to be estimated by ML)
#' @param X nxK model matrix of the convex combination
#' @param Z nxp model matrix
#' @param le nx1 vector of 1st derivatives of the log-likelihood with respect to eta
#' @param lee nx1 vector of 2nd derivatives of the log-likelihood with respect to eta
#' @param leee nx1 vector of 3rd derivatives of the log-likelihood with respect to eta
#' @param d1b (pK)xm matrix of derivatives of the Hessian with respect to the smoothing parameters
#' @param deriv 0: just grad and Hess, 1: first deriv of Hess
#'
#' @name convertDerivStackPositive
#' @rdname convertDerivStackPositive
#' @export
#' 
convertDerivStackPositive <- function(beta, X, Z,
                                      le, lt, 
                                      lee, let, ltt, 
                                      leee = 0, leet = 0, lett = 0, lttt = 0, 
                                      d1b = 0, deriv = 0) {
  
  d1H <- NULL ## default
  
  p <- ncol(Z); n <- nrow(X); K <- ncol(X); nth <- length(lt)
  beta <- matrix(beta, nrow = p, ncol = K)
  P <- p * K
  nu <- Z %*% beta
  a <- exp(nu)
  eta <- rowSums(a * X)
  
  ## the gradient...
  en <- X * a
  ln <- le * en
  lb <- c(as.vector(t(Z) %*% ln), lt)
  
  ## the Hessian...
  i2 <- mgcv::trind.generator(K)$i2
  lnn <- list()
  coun <- 0
  for (rr in 1:K) for (ss in rr:K) {
    coun <- coun + 1
    lnn[[coun]] <- lee * en[, rr] * en[, ss]
    if (rr == ss) lnn[[coun]] <- lnn[[coun]] + le * en[, rr] # note enn = en
  }
  lnn <- do.call(cbind, lnn)
  lbb <- lapply(1:K, function(x) vector("list", K))
  for (rr in 1:K) for (ss in rr:K) {
    lbb[[rr]][[ss]] <- t(Z * lnn[, i2[rr, ss]]) %*% Z
    if (ss > rr) lbb[[ss]][[rr]] <- lbb[[rr]][[ss]]
  }
  lbb <- lapply(lbb, function (x) do.call(rbind, x))
  lbb <- matrix(do.call(cbind, lbb), nrow = P)
  ltb <- list()
  coun <- 0
  for (rr in 1:K) {
    coun <- coun + 1
    ltb[[coun]] <- t(let) %*% (en[, rr] * Z)
  }
  ltb <- do.call(cbind, ltb)
  lbb <- rbind(cbind(lbb, t(ltb)), cbind(ltb, ltt))
  
  if (deriv == 1) { ## only store diagonal of d1H
    
    stop("diagonal elements of d1H only not provided, only full d1H implemented")
    
  } ## if deriv == 1
  
  if (deriv == 2) { ## store full d1H
    
    lnnn <- list()
    coun <- 0
    for (rr in 1:K) for (ss in rr:K) for (tt in ss:K) {
      coun <- coun + 1
      term1 <- leee * en[, rr] * en[, ss] * en[, tt]
      term2 <- term3 <- 0
      if (rr == ss) term2 <- term2 + en[, rr] * en[, tt]
      if (rr == tt) term2 <- term2 + en[, rr] * en[, ss]
      if (ss == tt) term2 <- term2 + en[, ss] * en[, rr]
      term2 <- lee * term2
      if (rr == ss & ss == tt) term3 <- le * en[, rr] # note ennn = enn = en
      lnnn[[coun]] <- term1 + term2 + term3
    }
    lnnn <- do.call(cbind, lnnn)
    
    m <- ncol(d1b)
    d1br <- d1b[1:(p * K), , drop = FALSE]
    d1tr <- d1b[(nrow(d1b) - nth + 1):nrow(d1b), , drop = FALSE]
    
    d1H <- list()
    i2 <- trind.generator(nth)$i2
    for (l in 1:m) {
      
      d1brl <- matrix(d1br[, l], nrow = p, ncol = K)
      
      i3 <- trind.generator(K)$i3
      lbbr <- lbbt <- lapply(1:K, function(x) vector("list", K))
      lbbt <- lapply(1:nth, function(x) lbbr)
      for (rr in 1:K) for (ss in rr:K) {
        
        A <- matrix(leet * en[, rr] * en[, ss], nrow = n)
        if (rr == ss) A <- A + let * en[, rr] 
        
        V1 <- rowSums((Z %*% d1brl) * lnnn[, i3[rr, ss, ]])
        V2 <- A %*% d1tr[, l, drop = FALSE]
        V <- as.numeric(V1 + V2)
        lbbr[[rr]][[ss]] <- t(Z) %*% (Z * V)
        if (ss > rr) lbbr[[ss]][[rr]] <- lbbr[[rr]][[ss]]
        
        for (uu in 1:nth) {
          lbbt[[uu]][[rr]][[ss]] <- t(Z) %*% (Z * A[, uu])
          if (ss > rr) lbbt[[uu]][[ss]][[rr]] <- lbbt[[uu]][[rr]][[ss]]
        }
      }
      lbbr <- lapply(lbbr, function (x) do.call(rbind, x))
      lbbr <- matrix(do.call(cbind, lbbr), nrow = P)
      
      lbtr <- list()
      coun <- 0
      for (uu in 1:nth) {
        coun <- coun + 1
        lbbt[[uu]] <- lapply(lbbt[[uu]], function (x) do.call(rbind, x))
        lbbt[[uu]] <- matrix(do.call(cbind, lbbt[[uu]]), nrow = P)
        lbtr[[uu]] <- lbbt[[uu]] %*% as.vector(d1brl) 
        lbtr[[uu]] <- lbtr[[uu]] + as.numeric(t(Z) %*% (
          as.numeric(lett[, i2[uu, 1:nth], drop = FALSE] %*% d1tr[, l]) * en))
        uu <- uu + 1
      }
      lbtr <- do.call(cbind, lbtr)
      
      lttr <- matrix(0, nrow = nth, ncol = nth)
      for (uu in 1:nth) for(vv in uu:nth) {
        lttr[uu, vv] <- sum(t(Z) %*% (en * lett[, i2[uu, vv]]) * d1brl) + 
          lttt[uu, vv, ] %*% d1tr[, l]
        if (vv > uu) lttr[vv, uu] <- lttr[uu, vv]
      }
      
      d1H[[l]] <- rbind(cbind(lbbr, lbtr), cbind(t(lbtr), lttr))
      
    } ## for (l in 1:m)
  } ## if deriv == 2
  
  list(lb=lb,lbb=lbb,d1H=d1H)
}
