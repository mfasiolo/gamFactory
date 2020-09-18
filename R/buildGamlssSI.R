#'
#' Function for building new GAMLSS families
#' 
#' @description XXX.
#' @param np XXX.
#' @name buildGamlssSI
#' @rdname buildGamlssSI
#' @export buildGamlssSI
#'
buildGamlssSI <- function(fam, effInfo){
  
  availableDeriv <- min(fam$availableDeriv, 3)
  cdf <- fam$cdf
  rd <- fam$rd
  qf <- fam$qf
  residuals <- fam$residuals
  nam <- paste0(fam$nam, "SI")
  np <- fam$np
  postproc <- fam$postproc
  okLinks <- fam$links
  
  llkFam <- fam$llk
  
  checkExtra <- fam$checkExtra
  if( is.null(checkExtra) ){ checkExtra <- function(.ex) NULL }
  
  initialize <- expression({
    if ( is.null(start) ) { start <- rep(1, ncol(x)) }
   # if ( is.null(start) ) { start <- family$initFun(y = y, nobs = nobs, E = E, x = x, family = family) }
  }) 
  
  defLinks <- lapply(okLinks, "[[", 1) # Default link function(s)
  
  outFam <- function(link = defLinks, extra = fam$extra){
    
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
      
      derLev <- switch(as.character(deriv), "0" = 0, "1" = 2, "2" = 3, "3" = 3, "4" = 4)
      
      effType <- effInfo$type
      ne <- length( effType )
      
      eff <- list()
      pen <- list()
      kk <- 1
      for(ii in 1:ne){
        if( effType[ii] == "standard" ){
          eff[[ii]] <- buildStandardEffect( X[ , effInfo$iec[[ii]], drop = FALSE] )  
        } else {
          if( effType[ii] == "singleIndex" ){
            eff[[ii]] <- buildSingleIndexEffect(Xi = effInfo$extra[[ii]]$Xi, 
                                                splineDes = constrSplineDes(k = effInfo$extra[[ii]]$k, 
                                                                            m = effInfo$extra[[ii]]$m, 
                                                                            lim = effInfo$extra[[ii]]$xlim)) 
            pen[[kk]] <- pen_varSI(a = coef[effInfo$iec[[ii]][1:ncol(effInfo$extra[[ii]]$Xi)]], 
                                   x = effInfo$extra[[ii]]$Xi, 
                                   v = 1, 
                                   deriv = derLev)
            if(deriv > 1){
              pen[[kk]]$outer <- pen_varSI_outer(a = coef[effInfo$iec[[ii]][1:ncol(effInfo$extra[[ii]]$Xi)]], 
                                                 x = effInfo$extra[[ii]]$Xi, 
                                                 DaDr = d1b[effInfo$iec[[ii]][1:ncol(effInfo$extra[[ii]]$Xi)], , drop = FALSE])
            }
            pen[[kk]]$iec <- effInfo$iec[[ii]][1:ncol(effInfo$extra[[ii]]$Xi)]
          } else {
            stop("Don't know this effect type")
          }
        }
      }
      
      # Build linear predictor object
      olp <- buildMultiLP(eff = eff, iel = effInfo$iel, iec = effInfo$iec)
      olp <- olp$eval(param = coef, deriv = derLev)
      
      # Evaluate eta and mu
      etas <- olp$f
      mus <- lapply(1:np, function(.kk) family$linfo[[.kk]]$linkinv( etas[[.kk]] ))
      
      # Derivatives of llk w.r.t. mu
      DllkDmu <- llkFam(y = y, param = mus, deriv = derLev)
      
      ret <- list("l" = drop(crossprod(wt, DllkDmu$d0)))
      
      lamVar <- 100
      
      npen <- length(pen) 
      if( npen ){
        for(ii in 1:npen){
          ret$l <- ret$l - lamVar * pen[[ii]]$d0
        }
      }
      
      if( deriv ){
        
        # Derivatives of llk w.r.t. eta 
        DllkDeta <- DllkDMu_to_DllkDeta(DllkDMu = DllkDmu, etas = etas, mus = mus, family = family, wt = wt, deriv = derLev)
        
        # Derivatives of log-likelihood w.r.t. beta
        DllkDbeta <- der(olp, param = coef, llk = DllkDeta, deriv = min(derLev, 2))
        
        ret$lb <- DllkDbeta$d1
        ret$lbb <- DllkDbeta$d2
        
        if( npen ){
          for(ii in 1:npen){
            # print(str(pen))
            # print(coef[1:4])
            # print(summary(effInfo$extra[[1]]$Xi%*%coef[1:4]))
            ret$lb[pen[[ii]]$iec] <- ret$lb[pen[[ii]]$iec] - lamVar * pen[[ii]]$d1
            ret$lbb[pen[[ii]]$iec, pen[[ii]]$iec] <- ret$lbb[pen[[ii]]$iec, pen[[ii]]$iec] - lamVar * pen[[ii]]$d2
          }
        }
        
        # print("###")
        # print(ret$l)
        # print(range(etas[[1]]))
        # print(var(Xsi %*% coefSI))
        print(coef)
        
        # We want also derivatives w.r.t. smoothing parameters
        if( deriv > 1 ){
          
          ret$d1H <- DHessDrho(o = olp, llk = DllkDeta, DbDr = d1b)
          
          if( npen ){
            for(ii in 1:npen){
              for(kk in 1:length(ret$d1H)){
                ret$d1H[[kk]][pen[[ii]]$iec, pen[[ii]]$iec] <- ret$d1H[[kk]][pen[[ii]]$iec, pen[[ii]]$iec] - lamVar * pen[[ii]]$outer[[kk]]
              }
            }
          }
          
        }
        
      }
      
      return( ret )
      
    } ## end ll 
    
    structure(list(family = nam, 
                   ll = ll, 
                   link = paste(link), 
                   nlp = np,
                   tri = trind.generator( np ), ## symmetric indices for accessing derivative arrays
                   initialize = initialize, 
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


