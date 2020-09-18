#'
#' Function for building new GAMLSS families
#' 
#' @description XXX.
#' @param np XXX.
#' @name buildGAMLSS
#' @rdname buildGAMLSS
#' @export buildGAMLSS
#'
buildGAMLSS <- function(o){
  
  availableDeriv <- o$availableDeriv
  cdf <- o$cdf
  rd <- o$rd
  qf <- o$qf
  residuals <- o$residuals
  llk <- o$llk  
  nam <- o$nam
  np <- o$np
  postproc <- o$postproc
  okLinks <- o$links
  
  checkExtra <- o$checkExtra
  if( is.null(checkExtra) ){ checkExtra <- function(.ex) NULL }
  
  initFun <- o$initFun
  initialize <- expression({
    if ( is.null(start) ) { start <- family$initFun(y = y, nobs = nobs, E = E, x = x, family = family) }
  }) 
  
  defLinks <- lapply(okLinks, "[[", 1)

  outFam <- function(link = defLinks, extra = o$extra){
    
    # Saving extra parameters in .GlobalEnv environment
    assign(".extra", extra, envir = environment())
    
    getExtra <- function( ) get(".extra")
    putExtra <- function(extra) assign(".extra", extra, envir = environment(sys.function()))
    
    if( !is.null(checkExtra) ){ checkExtra( extra ) }
    
    # Preparing link function
    stats <- getStats(link = link, 
                      okLinks = okLinks, 
                      np = np, 
                      nam = nam)
    
    ll <- function(y, X, coef, wt, family, offset=NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
      ## function defining a gamlss model log lik. 
      ## deriv: 0 - eval
      ##        1 - grad and Hess
      ##        2 - diagonal of first deriv of Hess
      ##        3 - first deriv of Hess
      ##        4 - everything.
      lpi <- attr(X, "lpi") ## extract linear predictor index
      np <- length(lpi)
      n <- length(y)
      
      extra <- get(".extra")
      
      # If offset for a linear predictor is absent or  NULL, put it to zero
      for(jj in 1:np){
        if( jj > length(offset) || is.null(offset[[jj]]) ){
          offset[[jj]] <- 0
        }
      }

      etas <- lapply(1:np, function(.kk) drop(X[ , lpi[[.kk]], drop=FALSE] %*% coef[ lpi[[.kk]] ]) + offset[[.kk]])
      mus <- lapply(1:np, function(.kk) family$linfo[[.kk]]$linkinv( etas[[.kk]] ))
      
      llkDer <- llk(y = y, 
                    param = do.call("cbind", mus), 
                    deriv = switch(as.character(deriv), "0" = 0, "1" = 2, "2" = 3, "3" = 3, "4" = 4), 
                    extra = extra)  
      
      l0  <- drop(crossprod(wt, llkDer$d0(SUM = FALSE)))
      
      if (deriv > 0) {
        
        l1 <- wt  * do.call("cbind", llkDer$d1(SUM = FALSE))
        l2 <- wt  * do.call("cbind", llkDer$d2(SUM = FALSE))

        ig1 <- do.call("cbind", lapply(1:np, function(.kk) family$linfo[[.kk]]$mu.eta( etas[[.kk]])))
        g2  <- do.call("cbind", lapply(1:np, function(.kk) family$linfo[[.kk]]$d2link( mus[[.kk]])))
        
      }
      
      l3 <- l4 <- g3 <- g4 <- 0 ## defaults
      
      if (deriv>1) {
        ## the third derivatives
        ## order mmm mmr mmx mrr mrx mxx rrr rrx rxx xxx
        l3 <- wt  * do.call("cbind", llkDer$d3(SUM = FALSE)) 

        g3  <- do.call("cbind", lapply(1:np, function(.kk) family$linfo[[.kk]]$d3link( mus[[.kk]])))

      }
      
      if (deriv>3) {
        ## the fourth derivatives
        ## mmmm mmmr mmmx mmrr mmrx mmxx mrrr mrrx mrxx mxxx
        ## rrrr rrrx rrxx rxxx xxxx
        l4 <- wt  * do.call("cbind", llkDer$d4()) 
        
        g4  <- do.call("cbind", lapply(1:np, function(.kk) family$linfo[[.kk]]$d4link( mus[[.kk]])))
        
      }
      if (deriv) {
        i2 <- family$tri$i2 
        i3 <- family$tri$i3
        i4 <- family$tri$i4
        
        ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
        de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)
        
        ## get the gradient and Hessian...
        ret <- gamlss.gH(X,lpi,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                         d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
      } else { 
        ret <- list()
      }
      
      ret$l <- l0; 
      
      return( ret )
    } ## end ll 
    
    structure(list(family = nam, 
                   ll = ll, 
                   link = paste(link), 
                   nlp = np,
                   tri = trind.generator( np ), ## symmetric indices for accessing derivative arrays
                   initialize = initialize, 
                   initFun = initFun,
                   postproc = postproc, 
                   residuals = residuals,
                   qf = qf,
                   linfo = stats,
                   rd = rd, 
                   cdf = cdf,
                   d2link=1,  ## signals to fix.family.link that all done 
                   d3link=1,
                   d4link=1, 
                   putExtra = putExtra, 
                   getExtra = getExtra,
                   ls = 1, ## signals that ls not needed here
                   available.derivs = availableDeriv - 2,
                   discrete.ok = TRUE
    ), class = c("general.family","extended.family","family"))
  
  }
  
}


