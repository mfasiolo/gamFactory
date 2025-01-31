#'
#' Preparing link function(s) of a response family
#' 
#' @description XXX.
#' @param np XXX.
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
    if( grepl("logea", link[[ii]], fixed=TRUE) ) { lnam <- "logea" }
    if (lnam %in% oklinks[[ii]] || 
        (lnam == "logitab" && grepl("logitab", oklinks[[ii]], fixed=TRUE)) || 
        (lnam == "loginva" && grepl("loginva", oklinks[[ii]], fixed=TRUE)) || 
        (lnam == "logea" && grepl("logea", oklinks[[ii]], fixed=TRUE))     ){
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