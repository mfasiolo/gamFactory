#'
#' Function for setting up link functions
#' 
#' @description Extension of \code{mgcv::make.link} including additional links.
#' @param link a character string indicating the link function to be used.
#' @name makeLink
#' @rdname makeLink
#' @importFrom gsubfn strapplyc
#' @export makeLink
#'
makeLink <- function(link){
  
  # If link is in mgcv this will work ...
  out <- tryCatch(make.link(link), error = function(e) e)
  
  # Otherwise it is defined directly here
  if( inherits(out, "error") ){
    
    # Get limits (a, b) of logitab link
    if( grepl("logitab", link, fixed=TRUE) ){
      
      tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
      l1 <- tmp[1]
      l2 <- tmp[2]
      link <- "logitab"
      
    }
    
    switch(link, 
           logitab = {
             linkfun <- eval(parse(text=paste("function(mu) binomial()$linkfun((mu - ", 
                                              l1, ") / ", l2 - l1, ")", sep='')))
             mu.eta <- eval(parse(text=paste("function(eta) binomial()$mu.eta(eta) * ", l2 - l1, sep=''))) 
             linkinv <- eval(parse(text=paste("function(eta) binomial()$linkinv(eta) * ",
                                              l2 - l1, " + ", l1, sep='')))
             valideta <- function(eta) all( is.finite(eta) )
           }, 
           stop(gettextf("%s link not recognised", sQuote(link)), domain = NA))
    
    environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- 
      environment(valideta) <- asNamespace("stats")
    
    out <- structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
                          valideta = valideta, name = link), class = "link-glm")
    
  }
  
  return( out )
}