#'
#' Preparing link function(s) of a response family
#' 
#' @description Internal function used to check whether a link function is within
#'             the set of allowed link functions for a given response family. It called
#'              \code{fix_family_link} to ensure that the link function has the correct derivatives.
#' @param link A character vector of link function(s) to be used for the response family.
#' @param oklinks A list of character vectors, each containing the allowed link functions.
#' @param np Number of parameters of the family.
#' @param nam Name of the family.
#' @name prepare_link
#' @rdname prepare_link
#' @export prepare_link
#' @examples
#' library(gamFactory)
#' prepare_link(link = "logitab(1.2, 4)", "logitab(3, 5)", 1, "gaussian")
#' prepare_link(link = "loginva(1.2)", "loginva", 1, "gaussian")
#' 
prepare_link <- function(link, oklinks, np, nam){
  
  if (length(link) != np){ stop( paste0("family requires ", np, "links specified as character strings") ) }
  
  stats <- list()
  for (ii in 1:np) {
    lnam <- link[[ii]]
    if( grepl("logitab", link[[ii]], fixed=TRUE) ) { lnam <- "logitab" }
    if( grepl("loginva", link[[ii]], fixed=TRUE) ) { lnam <- "loginva" }
    # NOTE: "logea2" contains "logea" as a substring, so its check must come after (and
    # therefore override) the "logea" check below.
    if( grepl("logea", link[[ii]], fixed=TRUE) ) { lnam <- "logea" }
    if( grepl("logea2", link[[ii]], fixed=TRUE) ) { lnam <- "logea2" }
    if (lnam %in% oklinks[[ii]] ||
        any(lnam == "logitab" & grepl("logitab", oklinks[[ii]], fixed=TRUE)) ||
        any(lnam == "loginva" & grepl("loginva", oklinks[[ii]], fixed=TRUE)) ||
        any(lnam == "logea" & grepl("logea", oklinks[[ii]], fixed=TRUE))     ||
        any(lnam == "logea2" & grepl("logea2", oklinks[[ii]], fixed=TRUE))   ){
      stats[[ii]] <- make_link(link[[ii]])
    } else {
      stop(lnam, paste0(" link not available for ", ii, "-th parameter of ", nam, " family"))
    }
    fam <- structure(list(link = link[[ii]], canonical = "none", 
                          linkfun = stats[[ii]]$linkfun, mu.eta = stats[[ii]]$mu.eta), 
                     class = "family")
    fam <- fix_family_link(fam)
    stats[[ii]]$d2link <- fam$d2link
    stats[[ii]]$d3link <- fam$d3link
    stats[[ii]]$d4link <- fam$d4link
  }
  
  return(stats)
}