#'
#' Calculating derivatives of log-likelihood wrt regression coefficients
#'
#' @description Intended for use with stack effects. 
#' Converts derivatives of the log-likelihood with respect to
#' the non linear predictor eta into derivatives with respect to the
#' regression coefficients and smoothing parameters
#'
#' @param param regression coefficients
#' @param theta vector of other parameters of the response variable distribution (to be included in beta to be estimated by ML)
#' @param X nxK model matrix of the convex combination
#' @param Z nxp model matrix
#' @param le nx1 vector of 1st derivatives of the log-likelihood with respect to eta
#' @param lee nx1 vector of 2nd derivatives of the log-likelihood with respect to eta
#' @param leee nx1 vector of 3rd derivatives of the log-likelihood with respect to eta
#' @param d1b (p(K-1))xm matrix of derivatives of the Hessian with respect to the smoothing parameters (at the moment only one smoothing parameter)
#' @param deriv 0: just grad and Hess, 1: first deriv of Hess
#'
#' @name convertDerivStackPositive
#' @rdname convertDerivStackPositive
#' @export
#' 
convertDerivStackPositive <- function(param, theta, X, Z, le, lee, leee = 0, 
                              d1b = 0, deriv = 0) {
  
  d1H <- NULL ## default
  
  p <- ncol(Z); n <- nrow(X); K <- ncol(X)
  param <- matrix(param, nrow = p, ncol = K)
  P <- p * K
  nu <- Z %*% param
  a <- exp(nu)
  eta <- rowSums(a * X)
  
  ## the gradient...
  en <- X * exp(nu)
  ln <- le * en
  lb <- as.vector(t(Z) %*% ln)
  
  ## the Hessian...
  enn <- list()
  coun <- 1
  for (jj in 1:K) for (kk in jj:K) {
    if (jj == kk) {
      enn[[coun]] <- exp(nu[, jj]) * X[, jj] 
    } else enn[[coun]] <- rep(0, n)
    coun <- coun + 1
  }
  enn <- do.call(cbind, enn)
  i2 <- mgcv::trind.generator(K)$i2
  lnn <- list()
  coun <- 1
  for (rr in 1:K) for (ss in rr:K) {
    lnn[[coun]] <- lee * en[, rr] * en[, ss] + le * enn[, i2[rr, ss]]
    coun <- coun + 1
  }
  lnn <- do.call(cbind, lnn)
  lbb <- lapply(1:K, function(x) vector("list", K))
  for (rr in 1:K) for (ss in rr:K) {
    lbb[[rr]][[ss]] <- t(Z * lnn[, i2[rr, ss]]) %*% Z
    if (ss > rr) lbb[[ss]][[rr]] <- lbb[[rr]][[ss]]
  }
  lbb <- lapply(lbb, function (x) do.call(rbind, x))
  lbb <- matrix(do.call(cbind, lbb), nrow = P)
  
  if (deriv == 1) {
    ennn <- list()
    coun <- 1
    for (jj in 1:K) for (kk in jj:K) for (ll in kk:K) {
      if(jj == kk & kk == ll) {
        ennn[[coun]] <- exp(nu[, jj]) * X[, jj] 
      } else ennn[[coun]] <- rep(0, n)
      coun <- coun + 1
    }
    ennn <- do.call(cbind, ennn)
    
    i3 <- mgcv::trind.generator(K)$i3
    lnnn <- list()
    coun <- 1
    for (rr in 1:K) for (ss in rr:K) for (tt in ss:K) {
      lnnn[[coun]] <- leee * en[, rr] * en[, ss] * en[, tt] +
        lee * (enn[, i2[rr, tt]] * en[, ss] + 
                 enn[, i2[ss, tt]] * en[, rr] + 
                 enn[, i2[rr, ss]] * en[, tt]) +
        le * ennn[, i3[rr, ss, tt]]
      coun <- coun + 1
    }
    lnnn <- do.call(cbind, lnnn)
    
    d1H <- lapply(1:K, function(x) vector("list", K))
    for (rr in 1:K) for (ss in rr:K) {
      d1H[[rr]][[ss]] <- t(Z) %*% (Z * rowSums((Z %*% matrix(d1b, nrow = p, ncol = K)) * 
                                                 lnnn[, i3[rr, ss, 1:K]]))
      if (ss > rr) d1H[[ss]][[rr]] <- d1H[[rr]][[ss]]
    }
    d1H <- lapply(d1H, function (x) do.call(rbind, x))
    d1H <- matrix(do.call(cbind, d1H), nrow = P)
    d1H <- as.numeric(d1H)
    
    # The following calculates the third derivatives wrt regression coefficients
    # We do not implement it because it should not be required by gam.fit5
    # lbbb <- lapply(1:K, function(x) {
    #   lapply(1:K, function(y) vector("list", K))
    # })
    # for (rr in 1:K) for (ss in 1:K) for (tt in 1:K) {
    #   lbbb[[rr]][[ss]][[tt]] <- array(0, dim = c(p, p, p))
    #   for (jj in 1:p) for (kk in 1:p) for (ll in 1:p) {
    #     lbbb[[rr]][[ss]][[tt]][jj, kk, ll] <- 
    #       sum(lnnn[, i3[rr, ss, tt]] * Z[, jj] * Z[, kk] * Z[, ll])
    #   }
    # }
    # for (rr in 1:K) for (ss in 1:K) {
    #   lbbb[[rr]][[ss]] <- do.call(abind, list(lbbb[[rr]][[ss]], along = 3))
    # }
    # for (rr in 1:K) {
    #   lbbb[[rr]] <- do.call(abind, list(lbbb[[rr]], along = 2))
    # }
    # lbbb <- do.call(abind, list(lbbb, along = 1))
    # 
    # lbbbMat <- list()
    # coun <- 1
    # for (jj in 1:P) for (kk in jj:P) for (ll in kk:P) {
    #   lbbbMat[[coun]] <- lbbb[jj, kk, ll]
    #   coun <- coun + 1
    # }
    # l3 <- do.call(c, lbbbMat)
  }
  
  list(lb=lb,lbb=lbb,d1H=d1H)
}
