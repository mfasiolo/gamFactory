
##################
#'
#' Log-likelihood of a multinomial logistic model
#'
#' @description Derivatives of the multinomial log-likelihood w.r.t. the linear predictors
#'              (log-odds relative to a reference category 0). Ported from the derivative
#'              computation used internally by [mgcv::multinom], but re-expressed as
#'              derivatives w.r.t. \code{param} (which plays the role of \code{mu}) rather
#'              than the model coefficients, and returned in the list format used by
#'              [gamFactory::llk_gaussian] and friends.
#' @param y a vector of observed category labels, integers between 0 (the reference
#'          category) and \code{K = ncol(param)}.
#' @param param a matrix (or list) with one column (element) per non-reference category,
#'              containing the corresponding linear predictor (log-odds relative to category 0).
#' @param deriv integer between 0 and 4 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-4
#'              additionally return \code{d1}-\code{d4}.
#' @rdname llk_multinomial
#' @export llk_multinomial
#' @examples
#' library(gamFactory)
#' n <- 10; K <- 3
#' y <- sample(0:K, n, replace = TRUE)
#' param <- matrix(rnorm(n * K), n, K)
#'
#' # Derivatives of multinomial log-lik up to order 4
#' llk_multinomial(y = y, param = param, deriv = 4)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_multinomial(y = y, param = matrix(param, n, K, byrow = TRUE), deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_multinomial(y = y, param = matrix(param, n, K, byrow = TRUE), deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_multinomial(y = y, param = matrix(param, n, K, byrow = TRUE), deriv = 2)$d2))
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_multinomial(y = y, param = matrix(param, n, K, byrow = TRUE), deriv = 3)$d3))
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_multinomial(y = y, param = matrix(param, n, K, byrow = TRUE), deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = rnorm(K), ord = 1:4)
#' @export
#'
llk_multinomial <- function(y, param, deriv = 0, ...) {

  if (is.list(param)) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)

  K <- ncol(param)
  n <- length(y)

  if (nrow(param) == 1 && n > 1) {
    param <- matrix(param, n, K, byrow = TRUE)
  }

  y <- round(y)
  if (min(y) < 0 || max(y) > K) stop("y must be between 0 and ncol(param)")

  ee <- exp(param)
  beta <- 1 + rowSums(ee)
  alpha <- log(beta)

  # eta0 has K+1 columns: column 1 is the (fixed to zero) reference category
  eta0 <- cbind(0, param)
  d0 <- eta0[cbind(1:n, y + 1)] - alpha

  out <- list("d0" = d0)

  tri <- trind.generator(K)

  if (deriv > 0) {

    p <- ee / beta # category probabilities for categories 1..K

    d1 <- lapply(1:K, function(.kk) as.numeric(y == .kk) - p[ , .kk])
    out[["d1"]] <- d1

    if (deriv > 1) {

      b2 <- beta^2
      d2 <- list()
      ii <- 0
      for (i in 1:K) for (j in i:K) {
        ii <- ii + 1
        d2[[ii]] <- if (i == j) { -p[ , i] + ee[ , i]^2 / b2 } else { (ee[ , i] * ee[ , j]) / b2 }
      }
      out[["d2"]] <- d2

      if (deriv > 2) {

        b3 <- b2 * beta
        d3 <- list()
        ii <- 0
        for (i in 1:K) for (j in i:K) for (k in j:K) {
          ii <- ii + 1
          if (i == j && j == k) { # all same
            d3[[ii]] <- d2[[ tri$i2[i, i] ]] + 2 * ee[ , i]^2 / b2 - 2 * ee[ , i]^3 / b3
          } else if (i != j && j != k && i != k) { # all different
            d3[[ii]] <- -2 * (ee[ , i] * ee[ , j] * ee[ , k]) / b3
          } else { # two same, one different
            kk <- if (i == j) k else j
            d3[[ii]] <- d2[[ tri$i2[i, kk] ]] - 2 * (ee[ , i] * ee[ , j] * ee[ , k]) / b3
          }
        }
        out[["d3"]] <- d3

        if (deriv > 3) {

          b4 <- b3 * beta
          d4 <- list()
          ii <- 0
          for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
            ii <- ii + 1
            uni <- unique(c(i, j, k, l))
            nun <- length(uni)
            if (nun == 1) { # all equal
              d4[[ii]] <- d3[[ tri$i3[i, i, i] ]] + 4 * ee[ , i]^2 / b2 - 10 * ee[ , i]^3 / b3 + 6 * ee[ , i]^4 / b4
            } else if (nun == 4) { # all unequal
              d4[[ii]] <- 6 * ee[ , i] * ee[ , j] * ee[ , k] * ee[ , l] / b4
            } else if (nun == 3) { # 2 same 2 different
              d4[[ii]] <- d3[[ tri$i3[uni[1], uni[2], uni[3]] ]] + 6 * ee[ , i] * ee[ , j] * ee[ , k] * ee[ , l] / b4
            } else if (sum(uni[1] == c(i, j, k, l)) == 2) { # 2 unique (2 of each)
              d4[[ii]] <- d3[[ tri$i3[uni[1], uni[2], uni[2]] ]] - 2 * ee[ , uni[1]]^2 * ee[ , uni[2]] / b3 +
                6 * ee[ , i] * ee[ , j] * ee[ , k] * ee[ , l] / b4
            } else { # 3 of one, 1 of the other
              if (sum(uni[1] == c(i, j, k, l)) == 1) uni <- uni[2:1]
              d4[[ii]] <- d3[[ tri$i3[uni[1], uni[1], uni[2]] ]] - 4 * ee[ , uni[1]]^2 * ee[ , uni[2]] / b3 +
                6 * ee[ , i] * ee[ , j] * ee[ , k] * ee[ , l] / b4
            }
          }
          out[["d4"]] <- d4

        }

      }

    }

  }

  return( out )

}
