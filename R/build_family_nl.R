#' 
#' Function for building GAM families containing non-standard effects
#' 
#' @name build_family_nl
#' @rdname build_family_nl
#' @export build_family_nl
#' 
build_family_nl <- function(bundle, info, lamVar = 1e5, lamRidge = 1e-5){
  
  available_deriv <- min(bundle$available_deriv, 3)
  cdf <- bundle$cdf
  rd <- bundle$rd
  qf <- bundle$qf
  store <- bundle$store
  residuals <- bundle$residuals
  bundle_nam <- bundle$bundle_nam
  nam <- paste0(bundle$nam, "_nl")
  np <- bundle$np
  postproc <- bundle$postproc
  oklinks <- bundle$links
  
  llkFam <- bundle$llk
  
  check_extra <- bundle$check_extra
  if( is.null(check_extra) ){ check_extra <- function(.ex) NULL }
  
  initialize_bundle <- bundle$initialize

  initialize_internal <- function(y, nobs, E, x, family, offset){
    
    p <- ncol( x )
    unscaled <- attr(E,"use.unscaled")
    lpi <- attr(x, "lpi")
    
    si <- which( sapply(info$type, paste0, collapse = '') != "stand" )
    nsi <- length( si )

    nkk <- 1:ncol(x)
    start <- numeric( ncol(x) )
    if( nsi ){ # If there are nested effects we: 
      # a) Identify coefficients of inner vector (alpha)
      iec <- info$iec[si]
      dsi <- sapply(info$extra[si], function(.x) length(.x$si$alpha))
      kk <- do.call("c", lapply(1:nsi, function(.ii) iec[[.ii]][1:dsi[.ii]]))
      # b) set the corresponding elements of start to the values contained in info$extra
      alpha <- do.call("c", lapply(info$extra[si], function(.x) .x$si$alpha))
      start[ kk ] <- alpha
      # c) modify lpi so that family$initialize_bundle will initialize only the remaining coefficients
      #    (excluding the columns of x and E related to the inner coefficients alpha)
      lpi <- lapply(lpi, function(.x) .x[ !(.x %in% kk) ])
      lpi <- lapply(lpi, function(.x) sapply(.x, function(.x1) .x1 - sum(kk < .x1)))
      E <- E[ , -kk, drop = FALSE]
      x <- x[ , -kk, drop = FALSE]
      attr(x, "lpi") <- lpi
      nkk <- nkk[ -kk ]
    }
    
    start[nkk] <- family$initialize_bundle(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset, 
                                            jj = lpi, unscaled = unscaled)
    
    return( start )
    
  }
  
  initialize <- expression({
    if ( is.null(start) ) { 
      start <- family$initialize_internal(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset) 
    }
  }) 
  
  defLinks <- lapply(oklinks, "[[", 1) # Default link function(s)
  
  outFam <- function(link = defLinks, extra = bundle$extra){
    
    # Saving extra parameters in .GlobalEnv environment
    assign(".extra", extra, envir = environment())
    
    get_extra <- function( ) get(".extra")
    put_extra <- function(extra) assign(".extra", extra, envir = environment(sys.function()))
    
    if( !is.null(check_extra) ){ check_extra( extra ) }
    
    # Preparing link function
    stats <- prepare_link(link = link, 
                          oklinks = oklinks, 
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
      
      ne <- length( info$type )
      
      # Build list of effect and penalties
      eff <- .build_effects(X = X, info = info, outer = outDer)
      
      # Build linear predictors and evaluate them
      olp <- linpreds(eff = eff, iel = info$iel, iec = info$iec)
      olp <- olp$eval(param = coef, deriv = derLev)
      
      # Evaluate effect-specific (not smoothing) penalties and their derivatives
      pen <- .eval_penalties(eff = olp$eff, info = info, d1b = d1b, deriv = derLev, outer = outDer)
      pen_ridge <- .eval_ridge_penalties(eff = olp$eff, info = info, deriv = derLev)
      
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
          ret$l <- ret$l - lamVar * pen[[ii]]$d0 - lamRidge * pen_ridge[[ii]]$d0
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
            ret$lb[zz] <- ret$lb[zz] - lamVar * pen[[ii]]$d1 - lamRidge * pen_ridge[[ii]]$d1
            ret$lbb[zz, zz] <- ret$lbb[zz, zz] - lamVar * pen[[ii]]$d2 - lamRidge * pen_ridge[[ii]]$d2
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
    #   effType <- info$type
    #   ne <- length( effType )
    #   
    #   # Build list of effect and penalties
    #   tmp <- .buildEffects(X = X, coef = beta, info = info, d1b = NULL, deriv = 0, outer = FALSE)
    #   eff <- tmp$eff
    # 
    #   # Build linear predictors and evaluate them
    #   olp <- linpreds(eff = eff, iel = info$iel, iec = info$iec)
    #   olp <- olp$eval(param = coef, deriv = derLev)
    #   
    #   if (se) return(list(fit=s,se.fit=sef)) else return(list(fit=olp$f))
    # } ## predict
    
    structure(list(family = nam, 
                   bundle_nam = bundle_nam,
                   ll = ll, 
                   link = paste(link), 
                   linkinv = if(np == 1){ stats[[1]]$linkinv } else { NULL },
                   nlp = np,
                   tri = trind.generator( np ), ## symmetric indices for accessing derivative arrays
                   initialize = initialize, 
                   initialize_bundle = initialize_bundle, 
                   initialize_internal = initialize_internal,
                   postproc = postproc, 
                   residuals = residuals,
                   store = store,
                   #predict = predict,
                   qf = qf,
                   linfo = stats,
                   rd = rd, 
                   cdf = cdf,
                   d2link=1,  ## signals to fix.family.link that all done 
                   d3link=1,
                   d4link=1, 
                   put_extra = put_extra, 
                   get_extra = get_extra,
                   ls = 1, ## signals that ls not needed here
                   available.derivs = available_deriv - 2,
                   discrete.ok = FALSE
    ), class = c("general.family","extended.family","family"))
    
  }
  
  return(outFam)
  
}


