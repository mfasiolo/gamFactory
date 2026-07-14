#'
#' Simulate data for testing gam_nl() nested effects
#'
#' @description Analogous to \code{\link[mgcv]{gamSim}}, but simulating data for the
#'              non-standard ("nested") effects built by \code{\link{s_nest}}
#'              (\code{\link{trans_linear}}, \code{\link{trans_mgks}},
#'              \code{\link{trans_nexpsm}}) fitted via \code{\link{gam_nl}}. Returns a
#'              data frame containing the response, the covariates needed to fit the
#'              model, and the true individual effects used to simulate the response
#'              (so that a fit can be checked against the truth), following
#'              \code{\link[mgcv]{gamSim}}'s own convention of returning \code{f}, and
#'              per-effect \code{f0}, \code{f1}, ... columns alongside the data.
#' @param eg integer, \code{1}, \code{2} or \code{3}, selecting which example to
#'           simulate:
#'           \describe{
#'             \item{\code{1}}{Two independent single-index effects
#'                   (\code{\link{trans_linear}}), each with a quadratic outer shape,
#'                   plus two standard smooths and a cubic term.}
#'             \item{\code{2}}{An exponentially-smoothed effect
#'                   (\code{\link{trans_nexpsm}}) plus a single-index effect, each with a
#'                   quadratic outer shape, plus a cubic term.}
#'             \item{\code{3}}{A multivariate-kernel-smoothed effect
#'                   (\code{\link{trans_mgks}}) plus a single-index effect, plus a cubic
#'                   term.}
#'           }
#' @param n number of observations to simulate. Defaults to a value tuned to keep the
#'          corresponding \code{\link{gam_nl}} fit reasonably fast (\code{500} for
#'          \code{eg = 1}, \code{300} for \code{eg = 2}, \code{900} for \code{eg = 3}, the
#'          latter rounded up to the nearest perfect square since it is arranged on a
#'          regular 2D grid).
#' @param verbose if \code{TRUE}, print a short description of the simulated example and
#'                the \code{\link{gam_nl}} formula needed to fit it.
#' @return A data frame with:
#'         \itemize{
#'           \item \code{y}: the simulated response.
#'           \item plain covariates needed to fit the model (e.g. \code{x}, \code{x1}).
#'           \item matrix-valued covariates for the nested effect(s) (\code{fake},
#'                 \code{fake1}, in the format expected by \code{\link{s_nest}} - see its
#'                 documentation and the examples below).
#'           \item \code{f}: the true overall mean used to simulate \code{y} (before
#'                 adding noise).
#'           \item one \code{f_*} column per individual true effect (e.g. \code{f_si1},
#'                 \code{f_si2} for \code{eg = 1}; \code{f_nexp}, \code{f_si} for
#'                 \code{eg = 2}/\code{3}), useful for checking that a fitted term
#'                 recovers the corresponding truth (e.g. via
#'                 \code{predict(fit, type = "terms")}).
#'         }
#'         The true single-index direction(s) are not per-observation, so they are
#'         attached as an attribute instead: \code{attr(out, "alpha")} is a named list
#'         with one direction vector per single-index effect in the example (\code{alpha1},
#'         \code{alpha2} for \code{eg = 1}; \code{alpha} for \code{eg = 2}/\code{3}).
#' @name gamFactorySim
#' @rdname gamFactorySim
#' @export gamFactorySim
#' @examples
#' library(gamFactory)
#'
#' ## eg = 1: two single-index effects
#' set.seed(1)
#' dat <- gamFactorySim(eg = 1)
#' fit <- gam_nl(list(y ~ s(x) + s_nest(fake, trans = trans_linear(pord = 1)) +
#'                      s(x1) + s_nest(fake1, trans = trans_linear(pord = 1)) + ti(x, x1),
#'                    ~ s(x, k = 4)),
#'               data = dat, n_init = 200, n_eigen = 3)
#' plot(dat$f, predict(fit)[ , 1]); abline(0, 1, col = 2)
#'
#' ## eg = 2: exponential smoothing + single-index effect
#' set.seed(1)
#' dat <- gamFactorySim(eg = 2)
#' fit <- gam_nl(list(y ~ s(x) +
#'                      s_nest(fake1, trans = trans_nexpsm(S = diag(3))) +
#'                      s_nest(fake2, trans = trans_linear(pord = 1)),
#'                    ~ s(x)),
#'               data = dat)
#' plot(dat$f, predict(fit)[ , 1]); abline(0, 1, col = 2)
#'
#' ## eg = 3: multivariate kernel smooth + single-index effect
#' set.seed(1)
#' dat <- gamFactorySim(eg = 3)
#' fit <- gam_nl(list(y ~ s(x) +
#'                      s_nest(fake1, trans = trans_mgks()) +
#'                      s_nest(fake2, trans = trans_linear()),
#'                    ~ 1),
#'               data = dat)
#' plot(dat$f, predict(fit)[ , 1]); abline(0, 1, col = 2)
#'
gamFactorySim <- function(eg = 1, n = NULL, verbose = TRUE){

  eg <- as.integer(eg)
  if( !(eg %in% 1:3) ){ stop("eg must be 1, 2 or 3") }

  descr <- c(
    "eg = 1: two single-index effects (trans_linear), each with a quadratic outer\n  shape, plus two standard smooths and a cubic term.\n  Fit via: gam_nl(list(y ~ s(x) + s_nest(fake, trans = trans_linear(pord = 1)) +\n                        s(x1) + s_nest(fake1, trans = trans_linear(pord = 1)) + ti(x, x1),\n                        ~ s(x, k = 4)), data = dat)",
    "eg = 2: exponentially-smoothed effect (trans_nexpsm) plus a single-index effect,\n  each with a quadratic outer shape, plus a cubic term.\n  Fit via: gam_nl(list(y ~ s(x) + s_nest(fake1, trans = trans_nexpsm(S = diag(3))) +\n                        s_nest(fake2, trans = trans_linear(pord = 1)),\n                        ~ s(x)), data = dat)",
    "eg = 3: multivariate-kernel-smoothed effect (trans_mgks) plus a single-index\n  effect, plus a cubic term.\n  Fit via: gam_nl(list(y ~ s(x) + s_nest(fake1, trans = trans_mgks()) +\n                        s_nest(fake2, trans = trans_linear()),\n                        ~ 1), data = dat)"
  )
  if( verbose ){ message(descr[eg]) }

  out <- switch(eg,
                `1` = .gamFactorySim_eg1(n = if(is.null(n)) 500 else n),
                `2` = .gamFactorySim_eg2(n = if(is.null(n)) 300 else n),
                `3` = .gamFactorySim_eg3(n = if(is.null(n)) 900 else n))

  return( out )

}

########################
# eg = 1: two single-index effects, two standard smooths, one cubic term
.gamFactorySim_eg1 <- function(n, p = 3){

  Xsi <- matrix(rt(p * n, df = 95), n, p, byrow = TRUE)
  alpha1 <- (1:p) / sum(1:p)
  alpha1 <- alpha1 / sd(Xsi %*% alpha1)
  siv <- drop( Xsi %*% alpha1 )

  Xsi1 <- matrix(rt(p * n, df = 95), n, p, byrow = TRUE)
  alpha2 <- (1:p) / sum(1:p)
  alpha2 <- alpha2 / sd(Xsi1 %*% alpha2)
  siv1 <- drop( Xsi1 %*% alpha2 )

  x <- runif(n, -1, 1)
  x1 <- runif(n, -1, 1)

  f_si1 <- 1 + siv + 2 * siv^2
  f_si2 <- 1 + siv1 - 2 * siv1^2
  f_x <- x + x^2
  f_x1 <- -x1^3

  f <- 1 + f_si1 + f_si2 + f_x + f_x1
  y <- f + rnorm(n)

  out <- data.frame(y = y, x = x, x1 = x1, f = f,
                    f_si1 = f_si1, f_si2 = f_si2, f_x = f_x, f_x1 = f_x1)
  out$fake <- Xsi
  out$fake1 <- Xsi1

  attr(out, "alpha") <- list(alpha1 = alpha1, alpha2 = alpha2)

  return( out )

}

########################
# eg = 2: exponentially-smoothed effect + single-index effect + cubic term
.gamFactorySim_eg2 <- function(n, p = 3){

  xseq <- seq(-4, 4, length.out = n)
  z <- 10 + sin(xseq) + rnorm(n)

  Xi <- cbind(1, xseq, xseq^2)
  beta <- c(2, 0, -0.1)
  zsm <- expsmooth(y = z, Xi = Xi, beta = beta)$d0
  f_nexp <- zsm + 0.5 * zsm^2

  Xsi <- matrix(rnorm(p * n), n, p, byrow = TRUE)
  alpha <- (1:p) / sum(1:p)
  siv <- drop( Xsi %*% alpha )
  f_si <- 1 + siv + 2 * siv^2

  x <- runif(n, -1, 1)
  f_x <- -0.3 * x^3

  f <- f_nexp + f_si + f_x
  y <- f + rnorm(n, 0, 0.5)

  out <- data.frame(y = y, x = x, f = f, f_nexp = f_nexp, f_si = f_si, f_x = f_x, zsm = zsm)

  fake1 <- cbind(z, Xi)
  colnames(fake1) <- c("y", rep("x", ncol(Xi)))
  out$fake1 <- fake1
  out$fake2 <- Xsi

  attr(out, "alpha") <- list(alpha = alpha)

  return( out )

}

########################
# eg = 3: multivariate-kernel-smoothed effect + single-index effect + cubic term
.gamFactorySim_eg3 <- function(n, p = 3, n0 = 60){

  ngr <- round(sqrt(n))
  n <- ngr^2

  X0 <- cbind(runif(n0, -1, 1), runif(n0, -4, 4))
  trueF <- function(x) 3 * x[ , 1] + x[ , 2]^2
  z <- trueF(X0)

  xseq1 <- seq(-1, 1, length.out = ngr)
  xseq2 <- seq(-4, 4, length.out = ngr)
  X <- as.matrix(expand.grid(xseq1, xseq2))

  dist <- lapply(1:ncol(X), function(dd) t(sapply(1:nrow(X), function(ii) (X[ii, dd] - X0[ , dd])^2)))
  beta <- c(-log(1), -log(2))
  zsm <- mgks(y = z, dist = dist, beta = beta)$d0
  f_mgks <- 1 + 0.3 * zsm - 0.05 * zsm^2

  Xsi <- matrix(rnorm(p * n), n, p, byrow = TRUE)
  alpha <- (1:p) / sum(1:p)
  siv <- drop( Xsi %*% alpha )
  f_si <- 1 + siv + 2 * siv^2

  x <- runif(n, -1, 1)
  f_x <- -0.3 * x^3

  f <- f_mgks + f_si + f_x
  y <- f + rnorm(n, 0, 0.5)

  out <- data.frame(y = y, x = x, f = f, f_mgks = f_mgks, f_si = f_si, f_x = f_x, zsm = zsm)

  dist0 <- do.call("cbind", dist)
  fake1 <- cbind(matrix(z, n, nrow(X0), byrow = TRUE), dist0)
  colnames(fake1) <- c(rep("y", nrow(X0)), rep("d1", ncol(dist0) / 2), rep("d2", ncol(dist0) / 2))
  out$fake1 <- fake1
  out$fake2 <- Xsi

  attr(out, "alpha") <- list(alpha = alpha)

  return( out )

}
