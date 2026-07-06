#'
#' Bundle for multinomial logistic regression model
#'
#' @description Ported from [`mgcv::multinom`]. There are \code{K} linear predictors, one
#'              per non-reference category, each giving the log-odds relative to
#'              category \code{0} (the reference category, which has no linear predictor
#'              of its own).
#' @param K number of linear predictors, i.e. number of categories minus one.
#' @name bundle_multinomial
#' @rdname bundle_multinomial
#' @export
#'
bundle_multinomial <- function(K){
  force(K)

  out <- list(np = K,
              available_deriv = 4,
              llk = llk_multinomial,
              links = replicate(K, "identity", simplify = FALSE),
              nam = "multinomial",
              store = list("K" = K),
              bundle_nam = as.character(match.call()[[1]]),

              # Predicted category probabilities, ported from mgcv::multinom()$predict.
              # Needed both for type = "response" predictions (mgcv::predict.gam looks
              # for family$predict when nlp > 1) and internally by "residuals" below.
              predict = function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                                  beta = NULL, off = NULL, Vb = NULL) {

                if (is.null(eta)) {
                  lpi <- attr(X, "lpi")
                  if (is.null(lpi)) { lpi <- list(1:ncol(X)) }
                  K <- length(lpi)
                  nobs <- nrow(X)
                  eta <- matrix(0, nobs, K)
                  if (se) {
                    ve <- matrix(0, nobs, K)
                    ce <- matrix(0, nobs, K * (K - 1) / 2)
                  }
                  ii <- 0
                  for (i in 1:K) {
                    Xi <- X[ , lpi[[i]], drop = FALSE]
                    eta[ , i] <- Xi %*% beta[lpi[[i]]]
                    if (!is.null(off[[i]])) eta[ , i] <- eta[ , i] + off[[i]]
                    if (se) {
                      ve[ , i] <- drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], lpi[[i]]]) * Xi)))
                      if (i < K) for (j in (i + 1):K) {
                        ii <- ii + 1
                        ce[ , ii] <- drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], lpi[[j]]]) * X[ , lpi[[j]]])))
                      }
                    }
                  }
                } else {
                  se <- FALSE
                  K <- ncol(eta)
                }

                gamma <- cbind(1, exp(eta))
                tot <- rowSums(gamma)
                gamma <- gamma / tot # category probabilities

                if (se) {
                  vp <- gamma * 0
                  for (j in 1:(K + 1)) {
                    if (j == 1) { dp <- -gamma[ , -1, drop = FALSE] / tot } else {
                      dp <- -gamma[ , j] * gamma[ , -1, drop = FALSE]
                      dp[ , j - 1] <- gamma[ , j] * (1 - gamma[ , j])
                    }
                    vp[ , j] <- rowSums(dp^2 * ve)
                    ii <- 0
                    for (i in 1:K) if (i < K) for (k in (i + 1):K) {
                      ii <- ii + 1
                      vp[ , j] <- vp[ , j] + 2 * dp[ , i] * dp[ , k] * ce[ , ii]
                    }
                    vp[ , j] <- sqrt(pmax(0, vp[ , j]))
                  }
                  return(list(fit = gamma, se.fit = vp))
                }

                return( list(fit = gamma) )

              },

              residuals = function(object, type = c("deviance")) {
                type <- match.arg(type)
                y <- round(object$y)
                n <- length(y)
                p <- object$family$predict(object$family, eta = object$fitted.values)[[1]]
                pc <- apply(p, 1, function(.x) which(max(.x) == .x)[1]) - 1
                sgn <- rep(-1, n)
                sgn[pc == y] <- 1
                return( sgn * sqrt(-2 * log(pmax(.Machine$double.eps, p[cbind(1:n, y + 1)]))) )
              },

              rd = function(mu, wt, scale) {
                p <- exp(cbind(0, mu))
                p <- p / rowSums(p)
                cp <- t(apply(p, 1, cumsum))
                return( apply(cp, 1, function(.x) min(which(.x > runif(1)))) - 1 )
              },

              postproc = expression({
                multinom <- list()
                object$y <- round(object$y)
                multinom$nj <- tabulate(object$y + 1)
                multinom$n <- sum(multinom$nj)
                multinom$K <- length(multinom$nj) - 1
                multinom$gamma <- c(1, solve(diag(multinom$n / multinom$nj[-1], multinom$K) -
                                    matrix(1, multinom$K, multinom$K), rep(1, multinom$K)))
                multinom$gamma <- log(multinom$gamma / sum(multinom$gamma))
                object$null.deviance <- -2 * sum(multinom$gamma[object$y + 1])
              }),

              initialize = function(y, nobs, E, x, family, offset, jj, unscaled, weights){

                K <- family$store$K
                lpi <- attr(x, "lpi")
                use.unscaled <- !is.null(attr(E, "use.unscaled"))
                start <- rep(0, ncol(x))

                for (k in 1:K) {
                  yt1 <- 6 * as.numeric(y == k) - 3
                  x1 <- x[ , lpi[[k]], drop = FALSE]
                  e1 <- E[ , lpi[[k]], drop = FALSE] ## square root of total penalty
                  if (use.unscaled) {
                    qrx <- qr(rbind(x1, e1))
                    x1 <- rbind(x1, e1)
                    startji <- qr.coef(qr(x1), c(yt1, rep(0, nrow(E))))
                    startji[!is.finite(startji)] <- 0
                  } else {
                    startji <- penreg(x1, e1, yt1)
                  }
                  start[lpi[[k]]] <- startji
                }

                return( start )
              }
  )

  # Fixing the environment of all functions
  for(ii in 1:length(out)){
    if( class(out[[ii]]) == "function" ){
      environment(out[[ii]]) <- environment()
    }
  }

  return( out )
}
