#'
#' Building new GAM families for mgcv::gam
#' 
#' @description XXX.
#' @param np XXX.
#' @name build_family
#' @rdname build_family
#' @export build_family
#' 
#' fam <- build_family(bundle_gaussian())
#' library(MASS)
#' b <- gam(list(accel~s(times,k=20,bs="ad"),~s(times)),
#'          data=mcycle,family=fam)
#' summary(b) 
#' plot(b,pages=1,scale=0)
#'
build_family <- function(bundle){
  
  available_deriv <- bundle$available_deriv
  cdf <- bundle$cdf
  rd <- bundle$rd
  qf <- bundle$qf
  residuals <- bundle$residuals
  bundle_nam <- bundle$bundle_nam
  nam <- paste0(bundle$nam, "_lss")
  np <- bundle$np
  postproc <- bundle$postproc
  oklinks <- bundle$links
  
  llkFam <- bundle$llk
  
  check_extra <- bundle$check_extra
  if( is.null(check_extra) ){ check_extra <- function(.ex) NULL }
  
  initialize_internal <- bundle$initialize
  
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
      
      
      etas <- lapply(1:np, function(.kk) drop(X[ , lpi[[.kk]], drop=FALSE] %*% coef[ lpi[[.kk]] ]) + offset[[.kk]])
      mus <- lapply(1:np, function(.kk) family$linfo[[.kk]]$linkinv( etas[[.kk]] ))
      
      # Derivatives of llk w.r.t. mu
      DllkDmu <- llkFam(y = y, param = mus, deriv = derLev)
      
      ret <- list()
      if( deriv ){
        
        # Derivatives of llk w.r.t. etas (linear predictors) 
        DllkDeta <- DllkDMu_to_DllkDeta(DllkDMu = DllkDmu, etas = etas, mus = mus, family = family, wt = wt, deriv = derLev)
        
        # We want also derivatives w.r.t. beta and rho (smoothing parameters)
        ret <- gamlss.gH(X,lpi,DllkDeta$l1,DllkDeta$l2,family$tri$i2,l3=DllkDeta$l3,i3=family$tri$i3,l4=DllkDeta$l4,i4=family$tri$i4,
                           d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
  
      }
      ret$l0 <- wt * DllkDmu$d0
      ret$l <- sum(ret$l0) 
      
      return( ret )
      
    } ## end ll 
    
    structure(list(family = nam, 
                   bundle_nam = bundle_nam,
                   ll = ll, 
                   link = paste(link), 
                   nlp = np,
                   tri = trind.generator( np ), ## symmetric indices for accessing derivative arrays
                   initialize = initialize, 
                   initialize_internal = initialize_internal,
                   postproc = postproc, 
                   residuals = residuals,
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
