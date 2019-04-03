d3r <- function(r, deriv = 1) {
  
  # penalty matrix with smoothing parameter lam
  S <- diag(rep(1, P))
  lam <- exp(r)
  
  # Inner step: optimize pen l wrt beta
  beta0 <- rep(0, P)
  ob <- createSHwithStack()$initialize(beta0, theta, X, Z, y)
  
  # OPTIM
  optimization <- optim(rep(0, P), function(x) {
    - ob$derObj(x)$d0() + .5 * lam * t(x) %*% S %*% x
    }, gr = function(x) - ob$derObj(x)$d1() + lam * S %*% x, 
    method = "BFGS", hessian = TRUE)
  betaHat <- optimization$par
  invNegH <- solve(optimization$hessian)
  
  # Given bhat, calculate derivatives 
  
  # d1b is the first derivative of b wrt rho 
  d1b <- - lam * invNegH %*% S %*% betaHat
  
  # Evaluate at beta hat
  # QUESTO CODICE SI RIPETE! IN REALTA TI SERVE SOLO LA DERIVATA RISPETTO A NU, 
  # OPPURE LA DERIVATA TERZA RISPETTO A BETAHAT CHE PERO E INCASINATA PER IL DISCORSO DEGLI INDICI
  # PER ORA VEDIAMO SOLO SE FUNZIONA
  # POI QUESTA PARTE VA CANCELLATA
  
  nuHat <- apply(matrix(betaHat, nrow = p, ncol = K - 1), 2, function(b) Z %*% b)
  aHat <- apply(nuHat, 1, function(nuHat) c(1, exp( nuHat )) / (1 + sum(exp(nuHat)))) %>% t
  etaHat <- rowSums(aHat * X)
  muHat <- etaHat
  
  
  nuHat <- Z %*% matrix(betaHat, nrow = p)
  aHat <- cbind(1, exp(nuHat)) / (1 + rowSums(exp(nuHat)))
  etaHat <- rowSums(aHat * X)
  muHat <- etaHat

  pars <- cbind(muHat, 1, 1, 1)
  objSH <- createSH(y = y)$derObj(param = pars, deriv = 3) # l wrt eta
  
  # eta wrt nu
  objStack <- createStackEffect(X)$derObj(nuHat, deriv = 3)
  
  le <- objSH$d1(SUM = FALSE)[[1]]
  en <- objStack$d1()
  lee <- objSH$d2(SUM = FALSE)[[1]]
  enn <- objStack$d2()
  leee <- objSH$d3(SUM = FALSE)[[1]]
  ennn <- objStack$d3()
  
  i2 <- trind.generator(K - 1)$i2
  i3rst <- trind.generator(K - 1)$i3
  
  # First derivative of the Hessian wrt to the smoothing parameters
  p <- ncol(Z)
  n <- length(y)
  K <- ncol(X)
  param <- matrix(betaHat, nrow = p, ncol = K - 1)
  P <- length(param)
  nu <- Z %*% param
  a <- cbind(1, exp(nu)) / (1 + rowSums(exp(nu)))
  eta <- rowSums(a * X)
  
  # With identity link mu=eta
  mu <- eta
  tau <- theta[1]
  eps <- theta[2]
  phi <- theta[3]
  sig <- exp(tau)
  del <- exp(phi)
  
  pars <- cbind(mu, 1, 1, 1)
  
  objSH <- createSH(y = y)$derObj(param = pars, deriv = 3)  # l wrt eta
  objStack <- createStackEffect(X)$derObj(nu, deriv = 3)    # eta wrt nu
  
  l <- sum(objSH$d0())
  le <- objSH$d1(SUM = FALSE)[[1]]
  en <- objStack$d1()
  l1 <- as.vector(t(Z * le) %*% en)
  
  lee <- objSH$d2(SUM = FALSE)[[1]]
  enn <- objStack$d2()
  i2 <- trind.generator(K - 1)$i2
  lbb <- lapply(1:(K - 1), function(x) vector("list", K - 1))
  for (rr in 1:(K - 1)) {
    for (ss in rr:(K - 1)) {
      lbb[[rr]][[ss]] <-
        t(Z * (lee * en[, rr] * en[, ss] + le * enn[, i2[rr, ss]])) %*% Z
      if (ss > rr) {
        lbb[[ss]][[rr]] <- lbb[[rr]][[ss]]
      }
    }
  }
  lbb <- lapply(lbb, function (x) do.call(rbind, x))
  l2 <- matrix(do.call(cbind, lbb), nrow = P)
  
  leee <- objSH$d3(SUM = FALSE)[[1]]
  ennn <- objStack$d3()
  
  i3rst <- trind.generator(K - 1)$i3
  
  
  # Calculate lnnn
  lnnn <- list()
  coun <- 1
  for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) for (tt in ss:(K - 1)) {
    lnnn[[coun]] <- leee * en[, rr] * en[, ss] * en[, tt] +
      lee * (enn[, i2[rr, tt]] * en[, ss] + 
               enn[, i2[ss, tt]] * en[, rr] + 
               enn[, i2[rr, ss]] * en[, tt]) +
      le * ennn[, i3rst[rr, ss, tt]]
    coun <- coun + 1
  }
  lnnn <- do.call(cbind, lnnn)
  
  lbbr <- lapply(1:(K - 1), function(x) vector("list", K - 1))
  for (rr in 1:(K - 1)) {
    for (ss in rr:(K - 1)) {
      lbbr[[rr]][[ss]] <- t(Z) %*% (Z * rowSums((Z %*% matrix(d1b, nrow = p)) * 
                                                  lnnn[, i3rst[rr, ss, 1:(K - 1)]]))
      if (ss > rr) {
        lbbr[[ss]][[rr]] <- lbbr[[rr]][[ss]]
      }
    }
  }
  lbbr <- lapply(lbbr, function (x) do.call(rbind, x))
  lbbr <- matrix(do.call(cbind, lbbr), nrow = P)
  lbbr <- as.numeric(lbbr)
  
  d0 <- function() {
    return(as.vector(ob$derObj(betaHat)$d2()))
    
  }
  d1 <- function() {
    if(deriv < 1) { stop("deriv < 1") }
    return(lbbr)
  }
  
  return(list("d0" = d0, "d1" = d1))
}

