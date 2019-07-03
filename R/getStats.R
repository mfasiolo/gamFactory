#'
#' Preparing link function(s) of a response family
#' 
#' @description XXX.
#' @param np XXX.
#' @name getStats
#' @rdname getStats
#' @export getStats
#'
getStats <- function(link, okLinks, np, nam){
  
  if (length(link) != np){ stop( paste0("family requires ", np, "links specified as character strings") ) }
  
  stats <- list()
  for (ii in 1:np) {
    if (link[[ii]] %in% okLinks[[ii]]){
      stats[[ii]] <- makeLink(link[[ii]])
    } else {
      stop(link[[ii]], paste0(" link not available for ", ii, "-th parameter of ", nam, " family"))
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