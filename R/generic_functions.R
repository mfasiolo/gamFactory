#'
#' Generic derivative function
#'
#' @param o object we want to get the derivatives of.
#' @param ... arguments to be passed to methods.
#' @rdname der
#' @export der
der <- function(o, ...) UseMethod("der")