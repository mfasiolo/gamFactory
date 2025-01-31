#'
#' Function for setting up link functions
#' 
#' @description Extension of \code{mgcv::make.link} including additional links.
#' @param link a character string indicating the link function to be used.
#' @name make_link
#' @rdname make_link
#' @importFrom gsubfn strapplyc
#' @export make_link
#'
make_link <- function(link){
  
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
      
    # Get limit (a, \infinity) on 1/mu of loginva link
    if( grepl("loginva", link, fixed=TRUE) ){
        
        tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
        l1 <- tmp[1]
        link <- "loginva"
        
    }
    
    # Get limit (a, \infinity) on exp(mu) logea link
    if( grepl("logea", link, fixed=TRUE) ){
      
      tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
      l1 <- tmp[1]
      link <- "logea"
      
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
           loginva = {
             linkfun <-  eval(parse(text=paste("function(mu) log(1/mu -",l1,")")))
             mu.eta <- eval(parse(text=paste("function(eta) { ee <- exp(eta); -ee/(ee +",l1,")^2 }")))
             linkinv <- eval(parse(text=paste("function(eta) 1/(exp(eta) +",l1,")")))
             valideta <- function(eta) all( is.finite(eta) )
           }, 
           logea = {
             linkfun <-  eval(parse(text=paste("function(mu) log(exp(mu) - ",l1,")", sep='')))
             mu.eta <-  eval(parse(text=paste("function(eta) { ee <- exp(eta); ee/(ee +",l1,") }")))
             linkinv <- eval(parse(text=paste("function(eta) log(exp(eta) +",l1,")", sep='')))
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