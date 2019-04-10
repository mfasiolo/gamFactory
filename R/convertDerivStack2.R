#'
#' Calculating derivatives of log-likelihood wrt regression coefficients
#'
#' @description Intended for use with stack effects. 
#' Converts derivatives of the log-likelihood with respect to
#' the non linear predictor eta into derivatives with respect to the
#' regression coefficients and smoothing parameters
#'
#' @param X nxK model matrix of the convex combination
#' @param Z nxp model matrix
#' @param lEtaList List of vectors of log-likelihood derivatives with respect to eta
#'
#' @name convertDerivStack
#' @rdname convertDerivStack
#' @export
#' 
convertDerivStack2 <- function(X, Z, lEtaList) {
  
  derObj <- function(param, d1b = 0, deriv = 0) {
    
    d1H <- NULL ## default
    
    p <- ncol(Z); n <- nrow(X); K <- ncol(X)
    param <- matrix(param, nrow = p, ncol = K - 1)
    P <- p * (K - 1)
    nu <- Z %*% param
    a <- cbind(1, exp(nu)) / (1 + rowSums(exp(nu)))
    eta <- rowSums(a * X)
    
    if (deriv > 0) {
      ## the gradient...
      am <- a[, - 1, drop = FALSE]
      xm <- X[, - 1, drop = FALSE]
      en <- am * (xm - eta)
      le <- lEtaList$le
      ln <- le * en
      lb <- as.vector(t(Z) %*% ln)
    }
    
    if (deriv > 1) {
      ## the Hessian...
      lee <- lEtaList$lee
      enn <- list()
      coun <- 1
      for (jj in 1:(K - 1)) for (kk in jj:(K - 1)) {
        enn[[coun]] <- en[, jj] * (as.numeric(jj==kk) - am[, kk]) - 
          am[, jj] * en[, kk]
        coun <- coun + 1
      }
      enn <- do.call(cbind, enn)
      i2 <- mgcv::trind.generator(K - 1)$i2
      lnn <- list()
      coun <- 1
      for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
        lnn[[coun]] <- lee * en[, rr] * en[, ss] + le * enn[, i2[rr, ss]]
        coun <- coun + 1
      }
      lnn <- do.call(cbind, lnn)
      lbb <- lapply(1:(K - 1), function(x) vector("list", K - 1))
      for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
        lbb[[rr]][[ss]] <- t(Z * lnn[, i2[rr, ss]]) %*% Z
        if (ss > rr) lbb[[ss]][[rr]] <- lbb[[rr]][[ss]]
      }
      lbb <- lapply(lbb, function (x) do.call(rbind, x))
      lbb <- matrix(do.call(cbind, lbb), nrow = P)
    }
    
    if (deriv > 2) {
      # The first derivative of the negative Hessian wrt smoothing parameters
      leee <- lEtaList$leee
      ennn <- list()
      coun <- 1
      for (jj in 1:(K - 1)) for (kk in 1:(K - 1)) for (ll in 1:(K - 1)) {
        ennn[[coun]] <- enn[, i2[jj, ll]] * (as.numeric(kk == jj) - am[, kk]) - 
          en[, jj] * am[, kk] * (as.numeric(kk == ll) - am[, ll]) - 
          en[, kk] * am[, jj] * (as.numeric(jj == ll) - am[, ll]) - 
          am[, jj] * enn[, i2[kk, ll]]
        coun <- coun + 1
      }
      ennn <- do.call(cbind, ennn)
      
      i3 <- mgcv::trind.generator(K - 1)$i3
      lnnn <- list()
      coun <- 1
      for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) for (tt in ss:(K - 1)) {
        lnnn[[coun]] <- leee * en[, rr] * en[, ss] * en[, tt] +
          lee * (enn[, i2[rr, tt]] * en[, ss] + 
                   enn[, i2[ss, tt]] * en[, rr] + 
                   enn[, i2[rr, ss]] * en[, tt]) +
          le * ennn[, i3[rr, ss, tt]]
        coun <- coun + 1
      }
      lnnn <- do.call(cbind, lnnn)
      
      l1H <- lapply(1:(K - 1), function(x) vector("list", K - 1))
      for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
        l1H[[rr]][[ss]] <- t(Z) %*% (Z * rowSums((Z %*% matrix(d1b, nrow = p, ncol = K - 1)) * 
                                                   lnnn[, i3[rr, ss, 1:(K - 1)]]))
        if (ss > rr) l1H[[ss]][[rr]] <- l1H[[rr]][[ss]]
      }
      l1H <- lapply(l1H, function (x) do.call(rbind, x))
      l1H <- matrix(do.call(cbind, l1H), nrow = P)
      l1H <- as.numeric(l1H)
      l1H <- - l1H # Since it is the derivative of the negative hessian
      
      # The following calculates the third derivatives wrt regression coefficients
      # We do not implement it because it should not be required by gam.fit5
      # lbbb <- lapply(1:(K - 1), function(x) {
      #   lapply(1:(K - 1), function(y) vector("list", K -1))
      # })
      # for (rr in 1:(K - 1)) for (ss in 1:(K - 1)) for (tt in 1:(K - 1)) {
      #   lbbb[[rr]][[ss]][[tt]] <- array(0, dim = c(p, p, p))
      #   for (jj in 1:p) for (kk in 1:p) for (ll in 1:p) {
      #     lbbb[[rr]][[ss]][[tt]][jj, kk, ll] <- 
      #       sum(lnnn[, i3[rr, ss, tt]] * Z[, jj] * Z[, kk] * Z[, ll])
      #   }
      # }
      # for (rr in 1:(K - 1)) for (ss in 1:(K - 1)) {
      #   lbbb[[rr]][[ss]] <- do.call(abind, list(lbbb[[rr]][[ss]], along = 3))
      # }
      # for (rr in 1:(K - 1)) {
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
    
    # Provide accessors to derivatives
    d1 <- function() {
      if(deriv < 1) stop("deriv < 1")
      return(lb)
    }
    d2 <- function() {
      if(deriv < 2) stop("deriv < 2")
      return(lbb)
    }
    d1H <- function() {
      if(deriv < 3) stop("deriv < 3")
      return(l1H)
    }
    
    return(list("d1" = d1, "d2" = d2, "d1H" = d1H))
    
  } # derObj
  
  initialize <- function(X = NULL, Z = NULL, 
                         p = NULL, K = NULL, 
                         lEtaList) {
    n <- length(le)
    if (is.null(X)) X <- matrix(rnorm(n * K), nrow = n, ncol = K)
    if (is.null(Z)) Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
    return(convertDerivStack(X = X, Z = Z, lEtaList = lEtaList))
  }
  
  return(list("derObj" = derObj, "initialize" = initialize))
  
} # convertDerivStack
