#'
#' Settings for single-index initialisation
#'
#' @description Bundles the tuning parameters controlling the multi-start search used to
#'              initialise single-index (\code{trans_linear}) nested effects, so that they can be
#'              passed around as a single argument. Used by \link{gam_nl} and \link{build_family_nl}.
#'
#' @param n_init number of candidate directions tried when initialising each true
#'               single-index nested effect.
#' @param n_eigen number of eigenvector-based candidate directions (out of \code{n_init})
#'                used when initialising each single-index nested effect.
#' @param oversample oversampling factor used to spread the remaining candidate
#'                    directions over the sphere.
#' @param seed seed used to generate the candidate directions on the sphere, so that
#'             the search is reproducible. If \code{NULL}, R's ambient random number
#'             stream is used instead (so results depend on and affect the caller's
#'             random state).
#' @return A list, with class \code{si_init_control}, holding the arguments above.
#' @name si_init_control
#' @rdname si_init_control
#' @export si_init_control
#'
si_init_control <- function(n_init = 1000, n_eigen = 10, oversample = 10, seed = 1){

  out <- list("n_init" = n_init, "n_eigen" = n_eigen, "oversample" = oversample, "seed" = seed)

  class(out) <- "si_init_control"

  return( out )

}
