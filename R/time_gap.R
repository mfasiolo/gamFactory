#' Compute adjacent time gaps and fill the first gap with the mode
#'
#' @description
#' Computes the differences between consecutive observations of a time variable
#' within a data frame, and fills the first missing gap with the modal value of
#' the observed time gaps. The function supports both numeric and
#' \code{Date}/\code{POSIXct} time formats, with optional conversion to
#' \code{difftime} objects.
#'
#' @param dat A \code{data.frame} containing a time column.
#' @param time_col Character string giving the name of the time variable in \code{dat}.
#' @param new_col Character string specifying the name of the new column that
#'        will store the computed time gaps (default: \code{"time_gap"}).
#' @param units Character; the time units for the computed differences when
#'        \code{time_col} is of class \code{Date} or \code{POSIXct}.
#'        One of \code{"secs"}, \code{"mins"}, \code{"hours"},
#'        \code{"days"}, or \code{"weeks"} (default: \code{"secs"}).
#' @param return_difftime Logical; whether to return the result as a
#'        \code{difftime} vector (default: \code{FALSE}, returns numeric).
#'
#' @return
#' A data frame identical to \code{dat}, with an additional column named
#' \code{new_col} that contains the computed time gaps.  
#' The first gap is filled by the mode (most frequent value) of all valid
#' subsequent gaps. If there are no valid gaps, the first element is set to \code{NA}.
#'
#' @export

time_gap <- function(dat, time_col, new_col = "time_gap",
                     units = c("secs", "mins", "hours", "days", "weeks"),
                     return_difftime = FALSE) {
  units <- match.arg(units)
  stopifnot(is.data.frame(dat))
  stopifnot(time_col %in% names(dat))
  
  tvec <- dat[[time_col]]
  n <- length(tvec)
  if (n == 0L) stop("time column is empty.")
  if (n == 1L) {
    dat[[new_col]] <- NA_real_
    if (return_difftime) dat[[new_col]] <- structure(dat[[new_col]], class = "difftime", units = units)
    return(dat)
  }
  
  # Adjacent differences
  d_raw <- diff(tvec)
  if (inherits(d_raw, "difftime")) {
    gap_num <- as.numeric(d_raw, units = units)
  } else {
    gap_num <- as.numeric(d_raw)
  }
  gap_num <- c(NA_real_, gap_num)
  
  # ---- Simplified mode calculation (using table directly) ----
  x_ok <- gap_num[-1]
  x_ok <- x_ok[is.finite(x_ok)]
  if (length(x_ok) == 0L) {
    mode_val <- NA_real_
  } else {
    x_ok <- round(x_ok, digits = 0)  # e.g., round to whole seconds
    tab <- table(x_ok)
    # Ensure the smallest numeric value is chosen in case of ties:
    # explicitly sort by numeric key before which.max
    ord <- order(as.numeric(names(tab)))
    tab <- tab[ord]
    mode_val <- as.numeric(names(tab)[which.max(tab)])
  }
  gap_num[1] <- mode_val
  # -----------------------------------------------------------
  
  # Return type
  if (return_difftime) {
    gap_out <- structure(gap_num, class = "difftime", units = units)
  } else {
    gap_out <- gap_num
  }
  
  dat[[new_col]] <- gap_out
  dat
}
