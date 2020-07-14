#'
#' Build a linear predictor
#' 
#' @param iel vectors indicating to which linear predictor each component belongs to.
#' @param iec list of vectors indicating which regression coefficients belong to each component.
#'  
#' @rdname buildMultipleLP
#' @export buildMultipleLP
#'
buildMultipleLP <- function(o, iel, iec){
  
  # Number of linear predictors
  nlp <- length( unique(iel) )

  # List of vectors indicating which componets belong to each linear predictor
  ile <- lapply(1:nlp, function(ii) which(iel == ii))
  
  # Number of components
  nc <- length( o )
  
  # Evaluate all the linear predictors
  eval <- function( param ){
    
    # Looping over LPs, then components of that LP
    eta <- lapply(ile, 
                  function(ii){
                    Reduce("+", lapply(ii, function(kk) o[[kk]]$eval(param = param[ iec[[kk]] ])$f ))
                  })
    
    return( list("f" = eta) )
    
  }
  
  # Derivatives of the log-likelihood w.r.t. regression coefficients of all linear predictors
  derivLogLik <- function(param, llk, deriv = 1, outer = FALSE){
    
    # Extract the part of the log-likelihood relevant for (mixed) derivatives wrt linear predictors i1, i2, i3
    .subset_llk <- function(llk, i1, i2 = NULL, i3 = NULL, deriv = 1, index = NULL){
      if( !is.list(llk$d1) ){ return(llk) }
      if( deriv >= 1 ){
        llk$d1 <- llk$d1[[ i1 ]]
        if( deriv >= 2 ){
          llk$d2 <- llk$d2[[ index$i2[i1, i2] ]]
          if( deriv >= 3){
            llk$d3 <- llk$d3[[ index$i3[i1, i2, i3] ]]
          }
        }
      }
      return(llk)
    }
    
    indLP <- trind.generator(nlp)
    
    print(deriv)
    der <- lapply(1:nc, function(ii) o[[ii]]$derivLogLik(param = param[ iec[[ii]] ], 
                                                         llk = .subset_llk(llk = llk, 
                                                                           i1 = iel[[ii]], i2 = iel[[ii]], i3 = iel[[ii]], 
                                                                           deriv = deriv, index = indLP),
                                                         deriv = deriv, 
                                                         outer = ifelse(deriv > 1, TRUE, outer)))
    
    d1 <- do.call("c", lapply(der, "[[", "d1"))
    
    d2 <- NULL 
    if( deriv > 1){ 
      
      # Hessian blocks
      H <- vector(mode = "list", length = nc^2)
      for( ir in 1:nc ){
        for( ic in ir:nc ){
          Hess <- paste0(".Hess.", class(o[[ic]]), "_", class(o[[ir]]))
          H[[ic*(ir-1)*nc]] <- if(ir != ic){ 
            do.call(Hess, list("o1" = der[[ic]], "o2" = der[[ir]], 
                               "llk" = .subset_llk(llk = llk, i1 = iel[[ic]], i2 = iel[[ir]], deriv = 2, index = indLP)))
          } else { 
            der[[ir]]$d2 
          }
          if( ir != ic ){
            H[[ir+(ic-1)*nc]] <- t( H[[ic+(ir-1)*nc]] )
          }
        }
      } 
      
      # Build Hessian # INEFFICIENT
      for(ir in 1:nc){
        d2 <- rbind(d2, do.call("cbind", H[(1+(ir-1)*nc):(ir*nc)]))
      }
      
    }
    
    xtra <- if( outer ){ list("der" = der, "iel" = iel, "iec" = iec, "ile" = ile) } else { NULL }
    
    out <- list("d1" = d1, "d2" = d2, 
                "xtra" = xtra, 
                "type" = "multiLP")
    
    return( out )
    
  }
  
  out <- structure(list("eval" = eval, "derivLogLik" = derivLogLik, 
                        "nc" = nc, "iel" = iel, "iec" = iec, "ile" = ile), 
                   class = "multiLP")
  
  return( out )
  
}