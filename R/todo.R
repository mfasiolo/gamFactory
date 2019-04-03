library(gamFactory)
library(qgam)
library(dplyr)
library(abind)

set.seed(0)

# p <- ncol(Z)
# n <- length(y)
# K <- ncol(X)
# param <- matrix(param, nrow = p, ncol = K - 1)
# P <- length(param)
# nu <- Z %*% param
# a <- cbind(1, exp(nu)) / (1 + rowSums(exp(nu)))
# eta <- rowSums(a * X)
# mu <- eta # With identity link mu=eta
# pars <- cbind(mu, theta[1], theta[2], theta[3])


n <- 1e3
p <- 3
K <- 3
P <- p * (K - 1)
beta <- rnorm(P)
Z <- cbind(1, matrix(runif(n * p), nrow = n))
Z <- scale(Z, scale = FALSE)
Z <- Z[, - 1, drop = FALSE]
X <- matrix(runif(n * K), nrow = n)
nu <- Z %*% matrix(beta, nrow = p, ncol = K - 1)
a <- cbind(1, exp(nu)) / (1 + rowSums(exp(nu)))
eta <- rowSums(a * X)
mu <- eta
tau <- 1
eps <- 1
phi <- 1
sig <- exp(tau)
del <- exp(phi)
y <- mu + (del * sig) * sinh((1/del) * asinh(qnorm(runif(n))) + (eps/del))
theta <- c(tau, eps, phi)
S <- diag(rep(1, P))
r <- 1
lam <- exp(r)
ob <- createSHwithStack()$initialize(beta, theta, X, Z, y, S, lam)

der <- derivCheck(np = 10,
                  parSim = function(n){ matrix(rnorm(n*(p)), n, P) },
                  obj = ob,
                  ord = 1:3,
                  trans = function(.x) {
                    si <- sign(.x)
                    return( si * sqrt(abs(.x)) )
                  },
                  Z=Z,y=y,X=X,theta=theta,S=S,lam=lam)

ciao <- fdDeriv(derFunWrapper(d3r), param = - 3, ord = 1)$fd1
ciao
plot(ciao[, 1], ciao[, 2])
abline(0,1)
