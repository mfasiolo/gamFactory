#
# Function for building new GAMLSS families
# 
# @name buildGamlssSI
# @rdname buildGamlssSI
# @export buildGamlssSI
#
.buildGamlssSI <- function(fam, effInfo, lamVar = 100){
  
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
  
  initializeInner <- fam$initialize

  initFun <- function(y, nobs, E, x, family, offset){
    
    p <- ncol( x )
    unscaled <- attr(E,"use.unscaled")
    lpi <- attr(x, "lpi")
    
    si <- which( effInfo$type != "standard" )
    nsi <- length( si )

    start <- numeric( ncol(x) )
    if( nsi ){ # If there are nested effects we: 
      # a) Identify coefficients of inner vector (alpha)
      iec <- effInfo$iec[si]
      dsi <- sapply(effInfo$extra[si], function(.x) length(.x$si$alpha))
      kk <- do.call("c", lapply(1:nsi, function(.ii) iec[[.ii]][1:dsi[.ii]]))
      # b) set the corresponding elements of start to the values contained in effInfo$extra
      alpha <- do.call("c", lapply(effInfo$extra[si], function(.x) .x$si$alpha))
      start[ kk ] <- alpha
      # c) modify lpi so that family$initializeInner will initialize only the remaining coefficients
      #    (in the output initializeInner, the values corresponding to the alphas will be set to zero)
      lpi <- lapply(lpi, function(.x) .x[ !(.x %in% kk) ])
    }
    
    start <- start + family$initializeInner(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset, 
                                            jj = lpi, unscaled = unscaled)
    
    return( start )
    
  }
  
  initialize <- expression({
    if ( is.null(start) ) { start <- family$initFun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset) }
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
      outDer <- deriv > 1
      
      effType <- effInfo$type
      ne <- length( effType )
      
      # Build list of effect and penalties
      eff <- .build_effects(X = X, effInfo = effInfo, outer = outDer)
      
      # Build linear predictors and evaluate them
      olp <- linpreds(eff = eff, iel = effInfo$iel, iec = effInfo$iec)
      olp <- olp$eval(param = coef, deriv = derLev)
      
      # Evaluate effect-specifit (not ridge) penalties and their derivatives
      pen <- .eval_penalties(eff = olp$eff, effInfo = effInfo, d1b = d1b, deriv = derLev, outer = outDer)
      
      # Evaluate eta and mu
      etas <- olp$f
      mus <- lapply(1:np, function(.kk) family$linfo[[.kk]]$linkinv( etas[[.kk]] ))
      
      # Derivatives of llk w.r.t. mu
      DllkDmu <- llkFam(y = y, param = mus, deriv = derLev)
      
      ret <- list("l" = drop(crossprod(wt, DllkDmu$d0)))
      
      # Add penalties to likelihood
      npen <- length(pen) 
      if( npen ){
        for(ii in 1:npen){
          ret$l <- ret$l - lamVar * pen[[ii]]$d0
        }
      }
      
      if( deriv ){
        
        # Derivatives of llk w.r.t. etas (linear predictors) 
        DllkDeta <- DllkDMu_to_DllkDeta(DllkDMu = DllkDmu, etas = etas, mus = mus, family = family, wt = wt, deriv = derLev)
        
        # Derivatives of log-likelihood w.r.t. beta
        Dbeta <- DllkDbeta(olp, param = coef, llk = DllkDeta, deriv = min(derLev, 2))
        ret$lb <- Dbeta$d1
        ret$lbb <- Dbeta$d2
        
        # Add derivatives of penalties w.r.t. beta
        if( npen ){
          for(ii in 1:npen){
            zz <- pen[[ii]]$iec
            ret$lb[zz] <- ret$lb[zz] - lamVar * pen[[ii]]$d1
            ret$lbb[zz, zz] <- ret$lbb[zz, zz] - lamVar * pen[[ii]]$d2
          }
        }
        
        # We want also derivatives w.r.t. rho (smoothing parameters)
        if( outDer ){
          
          ret$d1H <- DHessDrho(o = olp, llk = DllkDeta, DbDr = d1b)
          
          # Add derivatives of penalties w.r.t. rho
          if( npen ){
            for(ii in 1:npen){
              for(kk in 1:length(ret$d1H)){
                zz <- pen[[ii]]$iec
                ret$d1H[[kk]][zz, zz] <- ret$d1H[[kk]][zz, zz] - lamVar * pen[[ii]]$outer[[kk]]
              }
            }
          }
          
        }
        
      }
      
      return( ret )
      
    } ## end ll 
    
    # predict <- function(family,se=FALSE,eta=NULL,y=NULL,
    #                     X=NULL,beta=NULL,off=NULL,Vb=NULL) {
    # 
    #   effType <- effInfo$type
    #   ne <- length( effType )
    #   
    #   # Build list of effect and penalties
    #   tmp <- .buildEffects(X = X, coef = beta, effInfo = effInfo, d1b = NULL, deriv = 0, outer = FALSE)
    #   eff <- tmp$eff
    # 
    #   # Build linear predictors and evaluate them
    #   olp <- linpreds(eff = eff, iel = effInfo$iel, iec = effInfo$iec)
    #   olp <- olp$eval(param = coef, deriv = derLev)
    #   
    #   if (se) return(list(fit=s,se.fit=sef)) else return(list(fit=olp$f))
    # } ## predict
    
    structure(list(family = nam, 
                   ll = ll, 
                   link = paste(link), 
                   nlp = np,
                   tri = trind.generator( np ), ## symmetric indices for accessing derivative arrays
                   initialize = initialize, 
                   initializeInner = initializeInner, 
                   initFun = initFun,
                   postproc = postproc, 
                   residuals = residuals,
                   #predict = predict,
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


