########################
# Multi-start backfitting initialization for single-index ("si") nested effects.
#
# For a gam_nl() model, each "special" (non-standard) effect has coefficients split into
# an inner block (info$iec[[idx]][1:na], the single-index direction "alpha") and an outer
# block (the outer spline coefficients "beta"). eff_si() (R/eff_si.R) computes the
# effective single-index value as ax = Xi %*% (alpha + a0) (or Xi %*% exp(alpha + a0) if
# positive_si), where a0 = info$extra[[idx]]$si$a0 is a fixed offset. So for any desired
# *effective* direction v, the coefficient values to place are simply v - a0. We exploit
# this to treat "resolved via search", "resolved via the historical default", and "the
# current candidate under test" uniformly (see .fit_reduced_si below).
#
# Only effects of type "si" (info$type[[idx]][1] == "si") with positive_si == FALSE get
# the sphere/eigenvector multi-start search; mgks, nexpsm and positive_si single-index
# effects keep the historical behaviour (use the construction-time default alpha as-is,
# leaving their already-constructed outer basis columns untouched).

# Population-sd helper, matching the convention used in smooth_construct_si_smooth_spec.R
# (si$alpha is rescaled there so that this equals 1 for the default direction).
.sd_n <- function(x) {
  n <- length(x)
  sd(x) * sqrt(n - 1) / sqrt(n)
}

# a0 (the fixed additive offset in alpha + a0) is a concept specific to "si" (trans_linear)
# effects (see smooth_construct_si_smooth_spec.R); it is never set for mgks/nexpsm effects
# (see smooth_construct_mgks_smooth_spec.R / smooth_construct_nexpsm_smooth_spec.R and
# eff_mgks.R / eff_nexpsm.R, none of which reference a0 at all). Treat a missing a0 as an
# additive zero, so that v = alpha + a0 / alpha = v - a0 both correctly degenerate to the
# raw alpha for those effect types.
.si_a0 <- function(extra) {
  a0 <- extra$si$a0
  if (is.null(a0)) { a0 <- rep(0, length(extra$si$alpha)) }
  a0
}

########################
# Build a matrix of n_init candidate directions for a single-index effect with inner
# design matrix Xi (n x r). Every row v is rescaled so that Var(Xi %*% v) == 1.
# Testing code
# Xi <- matrix(rnorm(1000 * 20), nrow = 1000, ncol = 20)
# Xi <- scale(Xi, scale = F)
# alpha <- rep(1, 20)
# 
# # Should give a warning
# v <- .si_search_candidates(Xi, alpha, n_init = 1000, n_eigen = min(ncol(Xi), 10), oversample = 10)
# 
# range(apply(v, 1, function(x) sd(Xi %*% x) * sqrt(999/1000)))
#
.si_search_candidates <- function(Xi, alpha, n_init = 1000, n_eigen = min(ncol(Xi), 10), oversample = 10) {

  d <- ncol(Xi)
  
  if(abs(.sd_n(Xi %*% alpha) - 1) > 1e-6) {
    warning("Initial alpha does not have unit variance; rescaling to unit variance.")
    alpha <- alpha / .sd_n(Xi %*% alpha)
  }
  
  if( any(abs(colMeans(Xi)) > 1e-6) ) {
   stop("Covariates must be centered (mean zero) for single-index search.")
  }

  # Initialisation provided by smooth_construct_si_smooth_spec.R. It is such that 
  # Var(Xi %*% alpha) == 1, so we can use it as-is as the first candidate.
  cand <- t(as.vector(alpha))
  
  if(n_init > 1){
    XtX <- crossprod(Xi)
    
    # Add eigenvectors of the covariates cross-product. These will be rescaled, but not rotated!
    n_eigen <- min(n_eigen, d)
    ev <- eigen(XtX, symmetric = TRUE)$vectors[ , 1:n_eigen, drop = FALSE]
    ev <- t(ev) / sqrt(colMeans( (Xi %*% ev)^2 ))
    
    cand <- rbind(cand, ev)
    
    # Add vectors on a half-sphere (to avoid the sign ambiguity of the single-index direction)
    # Makes sure that the variance of the single index is 1 for all candidates, i.e. Var(Xi %*% v) == 1.
    # This is done using the Cholesky factor of the covariates cross-product, so that the half-sphere is 
    # stretched to match the covariates' geometry. This rotates the simulated vectors, while we do not want
    # to rotate the eigenvectors. 
    n_left <- n_init - nrow(cand)
    if(n_left > 0) {
      vsim <- .grid_on_half_sphere(N = n_left, r = d, oversample = oversample)
      ch <- chol(crossprod(Xi))
      vsim <- t( backsolve(ch, t(vsim)) ) * sqrt(nrow(Xi))
      cand <- rbind(cand, vsim)
    }

  }

  return( cand )
}

########################
# Given the set of all special-effect indices (nested_idx), a named list `resolved`
# (as.character(idx) -> effective direction v) for the effects currently fixed / under
# test, and `reeval_idx` (the subset of `resolved` whose v may differ from the
# construction-time default and therefore need their outer basis recomputed - i.e. the
# "si" effects currently being searched, never mgks/nexpsm/positive_si), build the
# reduced x/E/lpi handed to family$initialize_bundle():
#  - effects in `reeval_idx` get their outer/basis columns overwritten via
#    basis$evalX(Xi %*% v)$X0, and their inner columns dropped;
#  - effects in `resolved` but not in `reeval_idx` (mgks, nexpsm, positive_si) keep their
#    already-constructed outer columns untouched (they were built at the
#    construction-time default, which is exactly what `resolved` holds for them), and
#    only their inner columns are dropped;
#  - effects not yet in `resolved` are dropped entirely (both inner and outer columns).
.build_reduced_design_for_si_init <- function(x, E, lpi, info, nested_idx, resolved, reeval_idx) {

  drop_idx <- integer(0)

  for (idx in nested_idx) {

    iec <- info$iec[[idx]]
    extra <- info$extra[[idx]]
    na <- length(extra$si$alpha)
    inner <- iec[1:na]
    outer <- iec[-(1:na)]

    key <- as.character(idx)
    if (is.null(resolved[[key]])) {
      drop_idx <- c(drop_idx, iec)
      next
    }

    drop_idx <- c(drop_idx, inner)

    if (idx %in% reeval_idx) {
      v <- resolved[[key]]
      ax <- drop( extra$si$X %*% v )
      x[ , outer] <- extra$basis$evalX(x = ax, deriv = 0)$X0
    }
  }

  keep_idx <- setdiff(seq_len(ncol(x)), drop_idx)

  lpi_new <- lapply(lpi, function(.x) match(.x[!(.x %in% drop_idx)], keep_idx))

  list(x = x[ , keep_idx, drop = FALSE],
       E = E[ , keep_idx, drop = FALSE],
       lpi = lpi_new,
       keep_idx = keep_idx)
}

########################
# Build the reduced design for the current `resolved` set, run family$initialize_bundle
# on it, and scatter the result into a full-length coefficient vector (zero everywhere
# by default, so any not-yet-resolved "si" effect contributes exactly zero to the linear
# predictor - see eff_si.R / linpreds.R).
.fit_reduced_si <- function(y, nobs, x, E, lpi, info, family, offset, weights, unscaled, nested_idx, resolved, reeval_idx) {

  reduced <- .build_reduced_design_for_si_init(x = x, E = E, lpi = lpi, info = info,
                                                nested_idx = nested_idx, resolved = resolved,
                                                reeval_idx = reeval_idx)

  attr(reduced$x, "lpi") <- reduced$lpi

  start_reduced <- family$initialize_bundle(y = y, nobs = nobs, E = reduced$E, x = reduced$x,
                                             family = family, offset = offset, jj = reduced$lpi,
                                             unscaled = unscaled, weights = weights)

  start <- numeric( ncol(x) )
  start[reduced$keep_idx] <- start_reduced

  # Need to subtract a0 from the resolved effective direction to get the actual alpha to place in the start vector.
  for (idx in nested_idx) {
    key <- as.character(idx)
    if (is.null(resolved[[key]])) next
    iec <- info$iec[[idx]]
    na <- length(info$extra[[idx]]$si$alpha)
    start[iec[1:na]] <- resolved[[key]] - .si_a0(info$extra[[idx]])
  }

  start
}

########################
# Top-level orchestrator: sequentially (one pass, no re-visits) resolve (i.e. find good initialitation for) 
# every true single-index ("si", non positive_si) effect via a ~n_init-candidate search, scored by
# the exact model log-likelihood after fitting the remaining coefficients via family$initialize_bundle(). 
# Effects of other nested types (mgks, nexpsm) and positive_si single-index effects use the construction-time initialisation, alpha, 
# provided by smooth_construct_*_smooth_spec.R.
# NOTE: code is looking for initialisation alpha* = alpha + a0, so the returned start vector has alpha = alpha* - a0 for each effect.
#       The helper function .si_a0() returns a0 = 0 for mgks/nexpsm/positive_si effects, so the returned start vector has 
#       alpha = alpha* for those effect types.
.backfit_si_initialize <- function(y, nobs, E, x, family, offset, weights, info,
                                    n_init = 1000, n_eigen = 10, oversample = 10) {

  lpi <- attr(x, "lpi")
  unscaled <- attr(E, "use.unscaled")

  is_special <- !sapply(info$type, function(.z) identical(.z[1], "stand"))
  nested_idx <- which(is_special)

  if (length(nested_idx) == 0) {
    return( family$initialize_bundle(y = y, nobs = nobs, E = E, x = x, family = family,
                                      offset = offset, jj = lpi, unscaled = unscaled, weights = weights) )
  }

  is_si <- sapply(nested_idx, function(idx) identical(info$type[[idx]][1], "si"))
  is_positive_si <- sapply(nested_idx, function(idx) isTRUE(info$extra[[idx]]$si$positive_si))

  si_search <- nested_idx[is_si & !is_positive_si]
  si_fixed <- nested_idx[!(is_si & !is_positive_si)]

  resolved <- list()
  for (idx in si_fixed) {
    resolved[[as.character(idx)]] <- info$extra[[idx]]$si$alpha + .si_a0(info$extra[[idx]])
  }

  for (idx in si_search) {

    Xi <- info$extra[[idx]]$si$X
    a0 <- info$extra[[idx]]$si$a0
    alpha <- info$extra[[idx]]$si$alpha
    cands <- .si_search_candidates(Xi = Xi, alpha = alpha + a0, n_init = n_init, n_eigen = n_eigen, oversample = oversample)

    best_score <- -Inf
    best_v <- cands[1, ]

    for (ii in seq_len(nrow(cands))) {

      resolved_try <- resolved
      resolved_try[[as.character(idx)]] <- cands[ii, ]

      start_try <- .fit_reduced_si(y = y, nobs = nobs, x = x, E = E, lpi = lpi, info = info,
                                    family = family, offset = offset, weights = weights,
                                    unscaled = unscaled, nested_idx = nested_idx, resolved = resolved_try,
                                    reeval_idx = si_search)

      score <- family$ll(y, x, start_try, weights, family, offset = offset, deriv = 0)$l

      if ( is.finite(score) && score > best_score ) {
        best_score <- score
        best_v <- cands[ii, ]
      }
    }

    resolved[[as.character(idx)]] <- best_v

  }

  .fit_reduced_si(y = y, nobs = nobs, x = x, E = E, lpi = lpi, info = info, family = family,
                   offset = offset, weights = weights, unscaled = unscaled,
                   nested_idx = nested_idx, resolved = resolved, reeval_idx = si_search)

}
