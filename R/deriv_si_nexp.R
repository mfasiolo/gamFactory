#' Linear transform followed by exponential smoothing and its derivatives (up to 3rd order)
#'
#' @description
#' Computes a nested structure composed of:
#' \itemize{
#'   \item an inner \strong{linear projection}: \eqn{z_t = X_{si,t}^\top \alpha_{si}}, and
#'   \item an outer \strong{adaptive exponential smoothing} process:
#'         \eqn{s_t = \omega_t s_{t-1} + (1 - \omega_t) z_t},
#'         where \eqn{\omega_t = \text{sigmoid}(X_{nexp,t}^\top \alpha_{nexp})}.
#' }
#' 
#' This function returns the smoothed values and their derivatives 
#' with respect to both parameter groups \eqn{(\alpha_{nexp}, \alpha_{si})}, 
#' up to the third order. It serves as a core computational component 
#' for nested single-index + exponential smoothing effects in \pkg{gamFactory}.
#'
#' @param X_si Numeric matrix of dimension \eqn{n \times p}, 
#'             containing the covariates for the inner linear transformation layer.
#' @param X_nexp Numeric matrix of dimension \eqn{n \times q}, 
#'               containing the covariates for the exponential smoothing layer.
#' @param param Numeric vector of length \eqn{p + q}, 
#'              storing parameters in the order \code{[alpha_nexp (q), alpha_si (p)]}.
#' @param times Optional integer vector specifying which time indices 
#'              to return results for (default: all time points).
#' @param deriv Integer (0–3); the derivative order to compute.
#'              \code{deriv = 0} returns only smoothed values,
#'              higher values include derivatives up to that order.
#' @param alpha_center Optional numeric vector of length \eqn{p} 
#'                     (used for centering inner linear transformation); 
#'                     set to \code{NULL} to disable centering.
#' @param Z0 Optional numeric scalar giving the initial smoothed state \eqn{s_0};
#'           if \code{NULL}, the recursion starts with \eqn{s_1 = z_1}.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\code{d0}}{Vector of smoothed values \eqn{s_t}.}
#'   \item{\code{d1}}{Matrix of first derivatives (\eqn{n \times (p+q)}).}
#'   \item{\code{d2}}{Matrix of second derivatives (\eqn{n \times ((p+q)(p+q+1)/2)}).}
#'   \item{\code{d3}}{Matrix of third derivatives (\eqn{n \times \text{choose}(p+q+2, 3)}).}
#' }
#'
#' @examples
#' # -------------------- Example: construct data --------------------
#' set.seed(1)
#' n <- 200
#' n_si <- 3   # p
#' n_nexp <- 2 # q
#' X_si   <- matrix(rnorm(n * n_si),   n, n_si)
#' X_nexp <- matrix(rnorm(n * n_nexp), n, n_nexp)
#' alpha_center <- rep(1, n_si)  # Optional centering for inner layer
#' 
#' # Parameter order: [alpha_nexp, alpha_si], length q + p
#' param <- rnorm(n_nexp + n_si)
#' 
#' # ==============================================================================
#' # Case A: positive_si = FALSE
#' # ==============================================================================
#' cat("\n>>> Running checks for standard logic (positive_si = FALSE)...\n")
#' 
#' obj_standard <- list(
#'   d0 = function(param) {
#'     sum(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                       param = param, alpha_center = alpha_center,
#'                       positive_si = FALSE, 
#'                       deriv = 0)$d0)
#'   },
#'   d1 = function(param) {
#'     tmp <- colSums(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                                  param = param, alpha_center = alpha_center,
#'                                  positive_si = FALSE,
#'                                  deriv = 1)$d1)
#'     as.list(tmp)
#'   },
#'   d2 = function(param) {
#'     tmp <- colSums(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                                  param = param, alpha_center = alpha_center,
#'                                  positive_si = FALSE,
#'                                  deriv = 2)$d2)
#'     as.list(tmp)
#'   },
#'   d3 = function(param) {
#'     tmp <- colSums(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                                  param = param, alpha_center = alpha_center,
#'                                  positive_si = FALSE,
#'                                  deriv = 3)$d3)
#'     as.list(tmp)
#'   }
#' )
#' 
#' # Numerical vs analytical derivative check
#' check_std <- gamFactory:::check_deriv(obj = obj_standard, param = param, ord = 1:3)
#' 
#' 
#' # ==============================================================================
#' # Case B: positive_si = TRUE
#' # ==============================================================================
#' cat("\n>>> Running checks for exponential mapping (positive_si = TRUE)...\n")
#' 
#' obj_positive <- list(
#'   d0 = function(param) {
#'     sum(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                       param = param, alpha_center = alpha_center,
#'                       positive_si = TRUE, 
#'                       deriv = 0)$d0)
#'   },
#'   d1 = function(param) {
#'     tmp <- colSums(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                                  param = param, alpha_center = alpha_center,
#'                                  positive_si = TRUE,
#'                                  deriv = 1)$d1)
#'     as.list(tmp)
#'   },
#'   d2 = function(param) {
#'     tmp <- colSums(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                                  param = param, alpha_center = alpha_center,
#'                                  positive_si = TRUE,
#'                                  deriv = 2)$d2)
#'     as.list(tmp)
#'   },
#'   d3 = function(param) {
#'     tmp <- colSums(deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                                  param = param, alpha_center = alpha_center,
#'                                  positive_si = TRUE,
#'                                  deriv = 3)$d3)
#'     as.list(tmp)
#'   }
#' )
#' 
#' # Numerical vs analytical derivative check
#' check_pos <- gamFactory:::check_deriv(obj = obj_positive, param = param, ord = 1:3)
#' 
#' 
#' # ==============================================================================
#' # -------------------- Consistency test for partial time evaluation --------------------
#' # ==============================================================================
#' cat("\n>>> Running consistency test for 'times' argument...\n")
#' time_seq <- sort(sample(1:n, 25))
#' 
#' xsm <- deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                      param = param,
#'                      alpha_center = alpha_center,
#'                      positive_si = TRUE,
#'                      deriv = 3)
#' 
#' xsm_2 <- deriv_si_nexp(X_si = X_si, X_nexp = X_nexp,
#'                        param = param,
#'                        times = time_seq,
#'                        alpha_center = alpha_center,
#'                        positive_si = TRUE,
#'                        deriv = 3)
#' 
#' stopifnot(identical(xsm$d0[time_seq], xsm_2$d0))
#' stopifnot(identical(xsm$d1[time_seq, ], xsm_2$d1))
#' stopifnot(identical(xsm$d2[time_seq, ], xsm_2$d2))
#' stopifnot(identical(xsm$d3[time_seq, ], xsm_2$d3))
#' 
#' cat(">>> All consistency tests passed!\n")
#'
#' @name deriv_si_nexp
#' @rdname deriv_si_nexp
#' @export deriv_si_nexp

deriv_si_nexp <- function(X_si, X_nexp, param, times = NULL, deriv = 0,
                          alpha_center = NULL, Z0 = NULL,
                          positive_si = FALSE) {
  
  # upper limit: 3rd derivative
  if (deriv > 3) { warning("deriv = ", deriv, " exceeds 3; truncating to 3"); deriv <- 3 }
  
  Z0_scalar <- if (is.null(Z0)) NA_real_ else as.numeric(Z0)
  
  out <- gamFactory:::.deriv_si_nexp_cpp(
    X_si  = X_si,
    X_nexp= X_nexp,
    param = param,
    deriv = as.integer(deriv),
    alpha_center = alpha_center,   
    Z0 = Z0_scalar,   # NA_real_
    positive_si  = as.logical(positive_si)
  )
  
  if (!is.null(times)) {
    out$d0 <- out$d0[times]
    if (deriv >= 1) out$d1 <- out$d1[times, , drop = FALSE]
    if (deriv >= 2) out$d2 <- out$d2[times, , drop = FALSE]
    if (deriv >= 3) out$d3 <- out$d3[times, , drop = FALSE]
  }
  
  if (deriv == 0) return(list(d0 = out$d0))
  if (deriv == 1) return(list(d0 = out$d0, d1 = out$d1))
  if (deriv == 2) return(list(d0 = out$d0, d1 = out$d1, d2 = out$d2))
  list(d0 = out$d0, d1 = out$d1, d2 = out$d2, d3 = out$d3)
}








# # ======================== OLD R-Version ========================
# deriv_si_nexp_R <- function(X_si, X_nexp, param, times = NULL, deriv = 0,
#                             alpha_center = NULL, Z0 = NULL,
#                             positive_si = FALSE) {
#   # ---- upper limit: 3rd derivative ----
#   if (deriv > 3) { warning("deriv = ", deriv, " exceeds 3; truncating to 3"); deriv <- 3 }
# 
#   n      <- nrow(X_si)
#   n_si   <- ncol(X_si)     # p
#   n_nexp <- ncol(X_nexp)   # q
#   m_par  <- n_nexp + n_si  # m = q + p
# 
#   sigmoid <- function(u) 1 / (1 + exp(-u))
#   mat <- function(x, nr, nc) { dim(x) <- c(nr, nc); x }
#   upper_tri_vec <- function(M){
#     v <- c(); for(i in seq_len(nrow(M))) for(j in i:ncol(M)) v <- c(v, M[i, j]); v
#   }
# 
#   # ---------------- Parameter order: [alpha_nexp (q), alpha_si (p)] ----------------
#   alpha_nexp <- param[1:n_nexp]
#   alpha_si   <- param[(n_nexp+1):(n_nexp+n_si)]
# 
#   if (positive_si) {
#     alpha_si <- exp(alpha_si)
#   }
# 
#   # Inner linear projection z_t
#   if(!is.null(alpha_center)) {
#     z <- drop(X_si %*% (alpha_si + alpha_center))
#   } else {
#     z <- drop(X_si %*% alpha_si)
#   }
# 
#   # Outer driver and weights
#   u  <- drop(X_nexp %*% alpha_nexp)
#   om <- sigmoid(u)
# 
#   # ----- containers -----
#   s <- numeric(n)
#   e <- matrix(0, n, n_si)              # ds/d alpha_si
#   d <- matrix(0, n, n_nexp)            # ds/d alpha_nexp
#   h <- array(0, c(n, n_si, n_nexp))    # d^2 s / (alpha_si, alpha_nexp)
#   g <- array(0, c(n, n_nexp, n_nexp))  # d^2 s / (alpha_nexp, alpha_nexp)
# 
#   if (deriv >= 3) {
#     T3 <- array(0, c(n, n_si,  n_nexp, n_nexp))
#     U3 <- array(0, c(n, n_nexp, n_nexp, n_nexp))
#   }
# 
#   # ----- initialization -----
#   if (is.null(Z0)) {
#     s[1] <- z[1]
#   } else {
#     s[1] <- om[1] * Z0 + (1 - om[1]) * z[1]
#   }
#   e[1,] <- X_si[1,]
# 
#   # ----- recursion -----
#   if (n >= 2) {
#     for (t in 2:n) {
#       xt  <- X_nexp[t, ]
#       xs  <- X_si[t, ]
#       w   <- om[t]
#       s1  <- w * (1 - w)
#       s2  <- s1 * (1 - 2 * w)
#       s3  <- s1 * (1 - 6 * w + 6 * w^2)
# 
#       s[t] <- w * s[t-1] + (1 - w) * z[t]
# 
#       e[t, ] <- w * e[t-1, ] + (1 - w) * xs
#       d[t, ] <- w * d[t-1, ] + (s[t-1] - z[t]) * s1 * xt
# 
#       h_prev <- mat(h[t-1, , ], n_si, n_nexp)
#       h_curr <- outer(e[t-1, ] - xs, s1 * xt)
#       h[t, , ] <- w * h_prev + h_curr
# 
#       g_prev <- mat(g[t-1, , ], n_nexp, n_nexp)
#       g[t, , ] <- w * g_prev +
#         s1 * ( outer(xt, d[t-1, ]) + outer(d[t-1, ], xt) ) +
#         (s[t-1] - z[t]) * s2 * outer(xt, xt)
# 
#       if (deriv >= 3) {
#         for (a in 1:n_nexp) {
#           for (b in 1:n_nexp) {
#             T3[t, , a, b] <- w * T3[t-1, , a, b] +
#               s1 * ( xt[a] * h_prev[, b] + xt[b] * h_prev[, a] ) +
#               (e[t-1, ] - xs) * (s2 * xt[a] * xt[b])
#           }
#         }
#         for (c in 1:n_nexp) {
#           U3[t, , , c] <- w * U3[t-1, , , c] +
#             (s1 * xt[c]) * g_prev +
#             s2 * xt[c] * ( outer(xt, d[t-1, ]) + outer(d[t-1, ], xt) ) +
#             s1 * ( outer(xt, g_prev[, c]) + outer(g_prev[, c], xt) ) +
#             d[t-1, c] * s2 * outer(xt, xt) +
#             (s[t-1] - z[t]) * s3 * xt[c] * outer(xt, xt)
#         }
#       }
#     }
#   }
# 
#   # ===================== CHAIN RULE POST-PROCESSING =====================
# 
#   # ---- combine d1 ----
#   e_final <- e
#   if (positive_si) {
#     e_final <- sweep(e, 2, alpha_si, "*")
#   }
#   d1_out <- if (deriv >= 1) cbind(d, e_final) else NULL
# 
#   # ---- combine d2 ----
#   d2_out <- NULL
#   if (deriv >= 2) {
#     Hlist <- vector("list", n)
#     for (t in 1:n) {
#       h_mat <- mat(h[t, , ],    n_si,   n_nexp)
#       g_mat <- mat(g[t, , ],    n_nexp, n_nexp)
# 
#       H_si_si <- matrix(0, n_si, n_si)
# 
#       if (positive_si) {
#         h_mat <- h_mat * alpha_si
#         diag(H_si_si) <- e_final[t, ]
#       }
# 
#       Hs_t  <- rbind(
#         cbind(g_mat,          t(h_mat)),
#         cbind(h_mat,          H_si_si)
#       )
#       Hlist[[t]] <- upper_tri_vec(Hs_t)
#     }
#     d2_out <- do.call(rbind, Hlist)
#   }
# 
#   # ---- combine d3 ----
#   d3_out <- NULL
#   if (deriv >= 3) {
#     ncol_d3 <- choose(m_par + 2, 3)
#     d3_out  <- matrix(0, nrow = n, ncol = ncol_d3)
# 
#     for (t in 1:n) {
#       col <- 1
#       for (j in 1:m_par) {
#         for (k in j:m_par) {
#           for (l in k:m_par) {
#             si_count <- sum(c(j, k, l) > n_nexp)
# 
#             if (si_count == 0) {
#               # all nexp
#               d3_out[t, col] <- U3[t, j, k, l]
#             } else if (si_count == 1) {
#               # one si
#               idx <- c(j, k, l)
#               si_pos <- which(idx > n_nexp)
#               si_idx <- idx[si_pos] - n_nexp
#               nexp_idx <- sort(idx[-si_pos])
#               a <- nexp_idx[1]; b <- nexp_idx[2]
# 
#               val <- T3[t, si_idx, a, b]
#               if (positive_si) val <- val * alpha_si[si_idx]
#               d3_out[t, col] <- val
# 
#             } else if (si_count == 2) {
#               if (positive_si) {
#                 idx <- c(j, k, l)
#                 si_pos <- which(idx > n_nexp)
#                 si_idx1 <- idx[si_pos[1]] - n_nexp
#                 si_idx2 <- idx[si_pos[2]] - n_nexp
#                 a <- idx[-si_pos]
# 
#                 if (si_idx1 == si_idx2) {
#                   d3_out[t, col] <- h[t, si_idx1, a] * alpha_si[si_idx1]
#                 }
#               }
#             } else if (si_count == 3) {
#               if (positive_si) {
#                 idx <- c(j, k, l)
#                 si_idx1 <- idx[1] - n_nexp
#                 si_idx2 <- idx[2] - n_nexp
#                 si_idx3 <- idx[3] - n_nexp
# 
#                 if (si_idx1 == si_idx2 && si_idx2 == si_idx3) {
#                   d3_out[t, col] <- e_final[t, si_idx1]
#                 }
#               }
#             }
#             col <- col + 1
#           }
#         }
#       }
#     }
#   }
# 
#   # ---- times ----
#   if (!is.null(times)) {
#     d0_out <- s[times]
#     if (deriv >= 1) d1_out <- d1_out[times, , drop = FALSE]
#     if (deriv >= 2) d2_out <- d2_out[times, , drop = FALSE]
#     if (deriv >= 3) d3_out <- d3_out[times, , drop = FALSE]
#   } else {
#     d0_out <- s
#   }
# 
#   # ---- return ----
#   if (deriv == 0) return(list(d0 = d0_out))
#   if (deriv == 1) return(list(d0 = d0_out, d1 = d1_out))
#   if (deriv == 2) return(list(d0 = d0_out, d1 = d1_out, d2 = d2_out))
#   return(list(d0 = d0_out, d1 = d1_out, d2 = d2_out, d3 = d3_out))
# }
# # ======================== Test R vs C++ ========================
# set.seed(1)
# n <- 200; p <- 3; q <- 2
# X_si   <- matrix(rnorm(n*p), n,p)
# X_nexp <- matrix(rnorm(n*q), n,q)
# param  <- rnorm(q+p)
# times  <- sample.int(n, 80)
# 
# # ==============================================================================
# # Casw 1: positive_si = FALSE
# # ==============================================================================
# cat(">>> Testing standard logic (positive_si = FALSE) ...\n")
# 
# ref_std <- deriv_si_nexp_R(X_si, X_nexp, param, times=NULL, deriv=3, alpha_center=NULL, Z0=NULL, positive_si=FALSE)
# ans_std <- deriv_si_nexp  (X_si, X_nexp, param, times=NULL, deriv=3, alpha_center=NULL, Z0=NULL, positive_si=FALSE)
# 
# stopifnot(max(abs(ref_std$d0 - ans_std$d0)) < 1e-8)
# stopifnot(max(abs(ref_std$d1 - ans_std$d1)) < 1e-8)
# stopifnot(max(abs(ref_std$d2 - ans_std$d2)) < 1e-8)
# stopifnot(max(abs(ref_std$d3 - ans_std$d3)) < 1e-8)
# 
# # check 'times'
# ans_sub_std <- deriv_si_nexp(X_si, X_nexp, param, times=times, deriv=3, positive_si=FALSE)
# stopifnot(all.equal(ref_std$d0[times], ans_sub_std$d0))
# 
# cat("✅ Passed!\n")
# 
# # ==============================================================================
# # Casw 2: positive_si = TRUE
# # ==============================================================================
# cat("\n>>> Testing exponential mapping (positive_si = TRUE) ...\n")
# 
# ref_pos <- deriv_si_nexp_R(X_si, X_nexp, param, times=NULL, deriv=3, alpha_center=NULL, Z0=NULL, positive_si=TRUE)
# ans_pos <- gamFactory:::deriv_si_nexp  (X_si, X_nexp, param, times=NULL, deriv=3, alpha_center=NULL, Z0=NULL, positive_si=TRUE)
# 
# stopifnot(max(abs(ref_pos$d0 - ans_pos$d0)) < 1e-8)
# stopifnot(max(abs(ref_pos$d1 - ans_pos$d1)) < 1e-8)
# stopifnot(max(abs(ref_pos$d2 - ans_pos$d2)) < 1e-8)
# stopifnot(max(abs(ref_pos$d3 - ans_pos$d3)) < 1e-8)
# 
# # check 'times'
# ans_sub_pos <- gamFactory:::deriv_si_nexp(X_si, X_nexp, param, times=times, deriv=3, positive_si=TRUE)
# stopifnot(all.equal(ref_pos$d0[times], ans_sub_pos$d0))
# 
# stopifnot(max(abs(ref_pos$d1[times,] - ans_sub_pos$d1)) < 1e-8)
# stopifnot(max(abs(ref_pos$d2[times,] - ans_sub_pos$d2)) < 1e-8)
# stopifnot(max(abs(ref_pos$d3[times,] - ans_sub_pos$d3)) < 1e-8)
# 
# cat("✅ Passed!\n")
# 
# # ==============================================================================
# # compare by visual check (first 5 rows)
# # ==============================================================================
# cat("\n--- Visual Check (positive_si = TRUE, first 5 rows) ---\n")
# cat("d2 (R vs C++):\n")
# print(cbind(R = ref_pos$d2[1:5, 1], Cpp = ans_pos$d2[1:5, 1]))
# 
# cat("\nd3 (R vs C++):\n")
# print(cbind(R = ref_pos$d3[1:5, 1], Cpp = ans_pos$d3[1:5, 1]))


