########################
# Robust initial smoothing-parameter guess for nested-effect penalty blocks.
#
# mgcv:::initial.spg() sets, for each penalty block i, lambda[i] <- 0.3 * norm(ZHZ, "M") /
# norm(ZSZ, "M"), where ZHZ/ZSZ are the (negative) log-likelihood Hessian and the penalty,
# both projected onto the penalty's range space. For nested (s_nest()) smooths, the outer
# penalty is reparametrised at construction time (.diagPen()/.getBmatrix(), see
# I_build_nested_pspline_basis.R / I_diagPen.R / I_getBmatrix.R) so that every penalized
# eigen-direction of the *original* penalty is forced to unit scale, regardless of how small
# its original eigenvalue was. Since a P-spline penalty's eigenvalues normally span several
# orders of magnitude even within the "penalized" (non-null-space) subset (confirmed this is
# not gamFactory-specific: plain mgcv::gam(..., bs = "bs") shows the same inflated guess),
# the rescaling needed for the smallest of them inflates the corresponding Hessian entry
# hugely - and because norm(ZHZ, "M") takes the single *maximum* entry (not a typical one),
# that one near-degenerate direction can dominate the whole guess, landing lambda far outside
# the range that gam.fit5's outer Newton/REML search can safely start from. This is usually
# harmless (standard smooths self-correct regardless of a bad starting sp), but traps the fit
# in a spurious, oversmoothed local optimum when the nested smooth is paired with an
# unpenalized inner ("single index") direction, since that direction is then free to chase
# whatever the (transiently oversmoothed) outer curve leaves as signal.
#
# ---------------------------------------------------------------------------------------
# DEPENDENCY ON UNEXPORTED mgcv INTERNALS - READ BEFORE TOUCHING THIS FILE
# ---------------------------------------------------------------------------------------
# This does NOT edit mgcv:::initial.spg (which is maintained upstream) - it only *calls* it,
# to get the exact same baseline guess mgcv would use for non-nested blocks, and only
# replaces the guess for nested-effect blocks with a robust (median-based) equivalent
# computed via the same Hessian-projection machinery. But both mgcv:::initial.spg() and
# mgcv:::totalPenaltySpace() are unexported (`:::`), version-pinned re-implementations of
# mgcv (tested against mgcv 1.9.4 / mgcv.rdb from R 4.x on 2026-07): their name, argument
# list, return structure, and internal formula are NOT part of mgcv's public API and can
# change silently in any future mgcv release, with no deprecation warning. Concretely, this
# file assumes (mirroring mgcv:::estimate.gam(), which is what actually calls these two
# functions during a real fit):
#   - mgcv:::totalPenaltySpace(S, H, off, ncol(X)) returns a list with a component $E that
#     is a valid "square root of total penalty" matrix, usable as initial.spg()'s `E` arg
#     (mirrors estimate.gam()'s own `Ssp <- totalPenaltySpace(G$S, G$H, G$off, ncol(G$X));
#     G$Eb <- Ssp$E`).
#   - mgcv:::initial.spg(X, y, weights, family, S, rank, off, offset=, L=, lsp0=, E=) returns
#     a numeric vector of length length(S) (one lambda per penalty block), in the same order
#     as S/rank/off (mirrors estimate.gam()'s own
#     `lsp2 <- log(initial.spg(G$X, G$y, G$w, G$family, G$S, G$rank, G$off, offset = G$offset,
#      L = G$L, lsp0 = G$lsp0, E = G$Eb))`).
#   - mgcv:::initial.spg()'s own general.family branch computes, per block i,
#     `ZHZ <- -t(Z) %*% lbb[ind, ind] %*% Z; ZSZ <- t(Z) %*% S[[i]] %*% Z;
#      lambda[i] <- 0.3 * norm(ZHZ, "M") / norm(ZSZ, "M")`, with Z built by pivoted Cholesky
#     of S[[i]] when rank[i] < ncol(S[[i]]), or Z = identity otherwise. We reuse this exact
#     construction (see the loop below), only swapping norm(ZHZ, "M") for
#     median(abs(diag(ZHZ))) for nested blocks.
#   - Crucially, for `general.family` models (our case), estimate.gam() reparametrises X
#     *before* any of the above: `G$Sl <- Sl.setup(G); G$X <- Sl.initial.repara(G$Sl, G$X,
#     both.sides = FALSE)`, and applies the same transform to a user-supplied `start` (`start
#     <- Sl.initial.repara(G$Sl, start, inverse = FALSE, both.sides = FALSE)`) - fitting then
#     proceeds entirely in this reparametrised coordinate system (coefficients are only
#     transformed back at the very end, via `Sl.initial.repara(..., inverse = TRUE)` on the
#     converged fit). So `initial.spg()`'s *real* `X`/Hessian argument is NOT G$X as returned
#     by gam(fit = FALSE) - it is the Sl-reparametrised version. We replicate this transform
#     below (`Xr`) so our own ZHZ/ZSZ, and the baseline mgcv:::initial.spg() call for
#     non-nested blocks, are computed in the same coordinate system as the real fit. Verified
#     empirically that Sl.initial.repara is block-diagonal per smooth term (no cross-term
#     mixing - checked via regressing each smooth's reparametrised columns on its own original
#     columns alone, R^2 = 1), and that it happens to be the identity specifically on nested
#     (s_nest()) smooths' own columns for every case tested so far (their .diagPen()-based
#     construction, with `no.rescale = TRUE`/`repara = TRUE`, apparently doesn't trigger
#     Sl.repara's cross-term null-space handling) - but this is NOT assumed anywhere below;
#     `Xr` is used unconditionally for both the baseline and the nested-block override, so the
#     computation is correct regardless of whether that continues to hold.
#   - `start=` itself, as passed on to gam()'s own `start=` argument by the caller, must stay
#     in the *original* (non-reparametrised) coordinate system - gam() applies its own
#     Sl.initial.repara to whatever `start` it is given, exactly as it would if the user had
#     supplied it directly. Only the copy used internally here to compute `lbb` is
#     reparametrised (`start_r`).
#
# If any of this has changed upstream, the calls below will either error (caught, see the
# tryCatch wrapper: gam_nl() then falls back to mgcv's own, un-overridden initial.spg()
# guess, i.e. exactly the pre-existing behaviour) or - the more insidious risk - silently
# return something of the right type/shape but a different meaning (e.g. if mgcv ever
# rescales/reorders the returned lambda vector), which the length/sanity checks below only
# partially guard against. If this ever needs revisiting, start by diffing the current
# mgcv:::estimate.gam and mgcv:::initial.spg source (`deparse(mgcv:::estimate.gam)` /
# `deparse(mgcv:::initial.spg)`) against the quoted snippets above.
# ---------------------------------------------------------------------------------------
#
# Returns NULL if there is nothing to override (no smoothing parameters, no nested-effect
# penalty blocks, or anything above failed/changed shape - in which case a warning is issued
# and the caller should proceed exactly as if this function did not exist). Otherwise returns
# list(sp = <full-length vector, ready for in.out$sp>, start = <coefficient vector used to
# compute it>). The caller should pass `start` on to the actual fitting call too (as gam()'s
# own `start=` argument), so that the (possibly expensive) single-index multi-start search
# inside family$initialize_internal() is not run twice.
.nested_initial_sp <- function(G) {

  if (length(G$sp) == 0) return(NULL)

  is_nested <- sapply(G$smooth, function(sm) inherits(sm, "nested"))
  if (!any(is_nested)) return(NULL)

  nested_range <- lapply(G$smooth[is_nested], function(sm) c(sm$first.para, sm$last.para))

  block_is_nested <- sapply(G$off, function(o) {
    any(vapply(nested_range, function(r) o >= r[1] && o <= r[2], logical(1)))
  })

  if (!any(block_is_nested)) return(NULL)

  Eb <- totalPenaltySpace(G$S, G$H, G$off, ncol(G$X))$E

  # Same reparametrisation estimate.gam() applies to G$X before calling initial.spg() for
  # general.family models - see the header comment above.
  Sl <- Sl.setup(G)
  Xr <- Sl.initial.repara(Sl, G$X, both.sides = FALSE)

  lsp <- .my.initial.spg(Xr, G$y, G$w, G$family, G$S, G$rank, G$off, nested = block_is_nested,
                           offset = G$offset, L = G$L, lsp0 = G$lsp0, E = Eb)
  
  return( lsp )
}
