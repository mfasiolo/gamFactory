#' Prepare data for single-index + exponential smoothing models (Chicago example)
#'
#' @description
#' This function will be used in \code{\link{chicago_data}()} in case of double nested effect,
#' as this type of effect requires a different structure of data.
#' It extracts the response variable, computes time gaps between consecutive 
#' observations, and constructs a combined predictor matrix that includes
#' the original covariates, an optional intercept, and the computed time-gap feature.
#'
#' @param dat A data frame containing the response \code{y}, covariate matrix \code{X},
#'            and a time variable. Typically the return value of \code{\link{chicago_data}()}.
#' @param time_col Character; name of the time column (default: \code{"time"}).
#' @param y_col Character; name of the response/target column (default: \code{"death"}).
#' @param gap_col Character; name of the output column storing computed time gaps
#'                (default: \code{"time_gap_num"}).
#' @param intercept Logical; whether to add an intercept column (a column of ones)
#'                  to the predictor matrix in expsmooth matrix(default: \code{FALSE}).
#' @param ... Additional arguments passed to \code{\link{time_gap}()} 
#'
#' @return
#' A data frame with two components:
#' \itemize{
#'   \item \code{y}: the response vector.
#'   \item \code{xw}: a matrix (class \code{AsIs}) combining predictors, 
#'         optional intercept, and time-gap column.
#' }
#'
#'
#' @export

pre_si_nexp_chicago <- function(dat,
                     time_col = "time",
                     y_col    = "death",
                     gap_col  = "time_gap_num",
                     intercept = TRUE,
                     ...) {
  stopifnot(is.data.frame(dat))
  if (!("X" %in% names(dat))) stop("dat$X is missing: please ensure the predictor matrix is stored in dat$X.")
  if (!(y_col %in% names(dat))) stop(sprintf("Response column '%s' does not exist.", y_col))
  if (!(time_col %in% names(dat))) stop(sprintf("Time column '%s' does not exist.", time_col))
  
  n <- nrow(dat)
  
  # 1) y
  y <- dat[[y_col]]
  
  # 2) compute time_gap
  dat_new <- time_gap(dat, time_col, new_col = gap_col, return_difftime = FALSE, ...)
  w <- dat_new[[gap_col]]
  
  # 3) X matrix (ensure it's a matrix)
  x <- dat$X
  if (!is.matrix(x)) x <- as.matrix(x)
  
  # 4) assemble XW: x + (optional intercept column) + time_gap
  if (isTRUE(intercept)) {
    XW <- cbind(x, intercept = rep(1, n), time_gap = w)
  } else {
    XW <- cbind(x, time_gap = w)
  }
  rownames(XW) <- rownames(dat)
  
  # doy
  # ⚠️ not sure if this is correct
  time_vec <- dat[[time_col]]
  doy <- (time_vec - 0.5) %% 365.25 
  
  # 5) return the combined structure
  out <- data.frame(
    death = y, 
    xw    = I(XW), 
    time  = time_vec, 
    doy   = doy 
  )
  out
}