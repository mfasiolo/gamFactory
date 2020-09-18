#'
#' Preparing response family function
#' 
#' @description Extension of \code{mgcv::fix.family.link} including additional links.
#' @param np XXX.
#' @name fixFamilyLink
#' @rdname fixFamilyLink
#' @export fixFamilyLink
#'
fixFamilyLink <- function(fam){
  
  if ( !inherits(fam,"family") ){ stop("fam not a family object") }
  if ( !is.null(fam$d2link) && !is.null(fam$d3link) && !is.null(fam$d4link) ){ return(fam) }
  
  # If fam$link is in mgcv this will work ...
  famO <- tryCatch(fix.family.link(fam), error = function(e) e)
  
  # Otherwise it is defined directly here
  if( inherits(famO, "error") ){
    
    link <- fam$link
    
    # Get limits (a, b) of logitab link
    if( grepl("logitab", link, fixed=TRUE) ){
      
      tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
      l1 <- tmp[1]
      l2 <- tmp[2]
      link <- "logitab"
      
    }
  
    switch(link, 
           logitab = { # Logit in (a, b)
             fam$d2link <- eval(parse(text=paste("function(mu) { mu <- (mu - ", 
                                      l1, ") / ", l2 - l1, "; (1/(1 - mu)^2 - 1/mu^2) / ", 
                                      (l2 - l1)^2, " }", sep='')))
             fam$d3link <- eval(parse(text=paste("function(mu) { mu <- (mu - ", 
                                                 l1, ") / ", l2 - l1, "; (2/(1 - mu)^3 + 2/mu^3) / ", 
                                                 (l2 - l1)^3, " }", sep='')))
             fam$d4link <- eval(parse(text=paste("function(mu) { mu <- (mu - ", 
                                                 l1, ") / ", l2 - l1, "; (6/(1 - mu)^4 - 6/mu^4) / ", 
                                                 (l2 - l1)^4, " }", sep='')))
           }, 
           stop(gettextf("%s link not recognised", sQuote(link)), domain = NA))
    
    environment(fam$d2link) <- environment(fam$d3link) <- environment(fam$d4link) <- environment(fam$linkfun)
    
  } else {
    
    return(famO)
    
  }
  
  return( fam )
  
}


