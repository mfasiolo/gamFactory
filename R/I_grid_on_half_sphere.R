########################
# Sample "well spaced out" points on the hyper-sphere of radius 1, in r dimension
# (actual sphere for r = 3).
# SEE ALSO testing code below.
# NOTE: we are actually sampling half-sphere because we flip the vectors so that
# the first element (first coordinate) is positive
# Algorithm:
# [1] Sample randomly M points on the sphere
# [2] Pick the first point, then one on that is the least correlated with it,
#     the one that is the least correlated with both, and so on, until you have
#     N < M points.
#
# `seed`: if not NULL, the random pool (step [1]) is generated deterministically from
# `seed` instead of from R's ambient random stream - the same (N, r, oversample, seed)
# always gives the same output. This does NOT touch R's own random state as seen by the
# caller: the current .Random.seed is saved before, and restored (or removed again, if it
# did not exist yet) after, so calling this with a fixed seed has no side effect on any
# other random draws happening before or after it in the same session.
.grid_on_half_sphere <- function(N, r = 3, oversample = 10, seed = NULL) {

  M <- oversample * N

  stopifnot(oversample >= 1)

  if( !is.null(seed) ){
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if( had_seed ){ old_seed <- get(".Random.seed", envir = .GlobalEnv) }
    on.exit({
      if( had_seed ){
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if( exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) ){
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }

  Q <- matrix(rnorm(M * r), M, r)
  Q <- Q / sqrt(rowSums(Q^2))
  
  keep <- integer(N)
  keep[1] <- 1L
  
  selected <- rep(FALSE, M)
  selected[1] <- TRUE
  
  # score[i] = similarity of i-th candidate 1st selected direction
  score <- abs(drop(Q %*% Q[1, ]))
  score[selected] <- Inf
  
  for (k in 2:N) {
    # Pick candidate with smallest maximum projective similarity
    next_id <- which.min(score)
    
    keep[k] <- next_id
    selected[next_id] <- TRUE
    
    # Only compare all candidates to the newly selected point
    new_sim <- abs(drop(Q %*% Q[next_id, ]))
    
    # Update similarity vector of all candidates with nearest selected point
    score <- pmax(score, new_sim)
    score[selected] <- Inf
  }
  
  out <- Q[keep, , drop = FALSE]
  
  # If first coordinate is negative, flip orientation of vector
  neg <- out[ , 1] < 0
  out[neg, ] <- -out[neg, ]
  
  out
}

#############################
#### TESTING CODE
#############################

# -------------------------------------------------------------------------
# 3D plot (AI generated)
# -------------------------------------------------------------------------
# library(plotly)
# .plot_projective_q_3d <- function(Q, show_antipodes = TRUE, show_lines = TRUE) {
#   stopifnot(ncol(Q) == 3)
#   
#   N <- nrow(Q)
#   
#   # Unit sphere mesh
#   theta <- seq(0, 2 * pi, length.out = 80)
#   phi <- seq(0, pi, length.out = 40)
#   
#   sx <- outer(cos(theta), sin(phi))
#   sy <- outer(sin(theta), sin(phi))
#   sz <- outer(rep(1, length(theta)), cos(phi))
#   
#   p <- plot_ly() |>
#     add_surface(
#       x = sx,
#       y = sy,
#       z = sz,
#       opacity = 0.12,
#       showscale = FALSE,
#       hoverinfo = "skip"
#     ) |>
#     add_markers(
#       x = Q[, 1],
#       y = Q[, 2],
#       z = Q[, 3],
#       text = paste0("q", seq_len(N)),
#       marker = list(size = 5),
#       name = "Projective directions"
#     )
#   
#   if (show_antipodes) {
#     p <- p |>
#       add_markers(
#         x = -Q[, 1],
#         y = -Q[, 2],
#         z = -Q[, 3],
#         text = paste0("-q", seq_len(N)),
#         marker = list(size = 3, opacity = 0.35),
#         name = "Antipodal equivalents"
#       )
#   }
#   
#   if (show_lines) {
#     for (i in seq_len(N)) {
#       p <- p |>
#         add_trace(
#           x = c(-Q[i, 1], Q[i, 1]),
#           y = c(-Q[i, 2], Q[i, 2]),
#           z = c(-Q[i, 3], Q[i, 3]),
#           type = "scatter3d",
#           mode = "lines",
#           line = list(width = 2),
#           opacity = 0.25,
#           showlegend = FALSE,
#           hoverinfo = "skip"
#         )
#     }
#   }
#   
#   p |>
#     layout(
#       scene = list(
#         xaxis = list(title = "q1", range = c(-1, 1)),
#         yaxis = list(title = "q2", range = c(-1, 1)),
#         zaxis = list(title = "q3", range = c(-1, 1)),
#         aspectmode = "cube"
#       ),
#       title = "Projective design on the unit sphere"
#     )
# }
# 
# ##### TEST
# Q <- .grid_on_half_sphere(N = 100, r = 3, oversample = 10)
# 
# .plot_projective_q_3d(Q)