#'
#' Build a linear predictor
#' 
#' @rdname buildLinearPredictor
#' @export buildLinearPredictor
#'
buildLinearPredictor <- function(o, KK){
  
  # Number of components
  nc <- length( o )
  
  eval <- function( param ){
    
    eta <- lapply(1:nc, function(ii) o[[ii]]$eval(param = param[ KK[[ii]] ])$f )
    eta <- Reduce("+", eta)
    
    return( list("f" = eta) )
    
  }
  
  derivLogLik <- function(param, llk, deriv = 1, outer = FALSE){
    
    der <- lapply(1:nc, function(ii) o[[ii]]$derivLogLik(param = param[ KK[[ii]] ], 
                                                         llk = llk, 
                                                         deriv = deriv, 
                                                         outer = ifelse(deriv > 1, TRUE, outer)))
    
    d1 <- do.call("c", lapply(der, "[[", "d1"))
    
    d2 <- NULL 
    if( deriv > 1){ 
    
      # Hessian blocks # INEFFICIENT should get only upper triangle
      jj <- 1
      H <- list()
      for( ir in 1:nc ){
        for( ic in 1:nc ){
          Hess <- paste0(".Hess.", class(o[[ic]]), "_", class(o[[ir]]))
          H[[jj]] <- if(ir != ic){ 
            do.call(Hess, list("o1" = der[[ic]], "o2" = der[[ir]], "llk" = llk))
          } else { 
            der[[ir]]$d2 
          }
          jj <- jj + 1
        }
      } 
      
      # Build Hessian # INEFFICIENT
      for(ir in 1:nc){
        d2 <- rbind(d2, do.call("cbind", H[(1+(ir-1)*nc):(ir*nc)]))
      }
      
    }
    
    xtra <- if( outer ){ list("der" = der, "KK" = KK) } else { NULL }
    
    out <- list("d1" = d1, "d2" = d2, 
                "xtra" = xtra, 
                "type" = "linearPredictor")
    
    return( out )
    
  }
  
  out <- structure(list("eval" = eval, "derivLogLik" = derivLogLik, 
                        "nc" = nc, "KK" = KK), 
                   class = "linearPredictor")
  
  return( out )
  
}
