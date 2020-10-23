#'
#' Preparing link function(s) of a response family
#' 
#' @description XXX.
#' @param np XXX.
#' @name getStats
#' @rdname getStats
#' @export getStats
#' @examples
#' library(gamFactory)
#' getStats(link = "logitab(1.2, 4)", "logitab(3, 5)", 1, "gaussian")
#' getStats(link = "loga(1.2)", "loga", 1, "gaussian")
#' 
getStats <- function(link, okLinks, np, nam){
  
  if (length(link) != np){ stop( paste0("family requires ", np, "links specified as character strings") ) }
  
  stats <- list()
  for (ii in 1:np) {
    lnam <- link[[ii]]
    if( grepl("logitab", link[[ii]], fixed=TRUE) ) { lnam <- "logitab" }
    if( grepl("loga", link[[ii]], fixed=TRUE) ) { lnam <- "loga" }
    if (lnam %in% okLinks[[ii]] || (lnam == "logitab" && grepl("logitab", okLinks[[ii]], fixed=TRUE)) || (lnam == "loga" && grepl("loga", okLinks[[ii]], fixed=TRUE))     ){
      stats[[ii]] <- makeLink(link[[ii]])
    } else {
      stop(lnam, paste0(" link not available for ", ii, "-th parameter of ", nam, " family"))
    }
    fam <- structure(list(link = link[[ii]], canonical = "none", 
                          linkfun = stats[[ii]]$linkfun, mu.eta = stats[[ii]]$mu.eta), 
                     class = "family")
    fam <- fixFamilyLink(fam)
    stats[[ii]]$d2link <- fam$d2link
    stats[[ii]]$d3link <- fam$d3link
    stats[[ii]]$d4link <- fam$d4link
  }
  
  return(stats)
}