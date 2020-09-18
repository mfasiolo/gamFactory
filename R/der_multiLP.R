#'
#' Derivatives of log-likelihood with multiple linear predictor
#' 
#' @rdname der.multiLP
#' @export der.multiLP
#' @export
der.multiLP <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv > 2 ) stop("Does not work with deriv > 2") 
  nlp <- o$nlp
  iec <- o$iec
  iel <- o$iel
  nc  <- o$nc
  
  if( is.null(param) ){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  }
  
  # Need to update the object
  if( is.null(o$param) || !identical(param, o$param) || o$deriv < deriv ){
    o <- o$eval(param = param, deriv = deriv)
  }
  
  subset_llk <- .create_subset_llk(llk, nlp) 
  
  der <- lapply(1:nc, function(ii) der(o$eff[[ii]], 
                                       param = param[ iec[[ii]] ], 
                                       llk = subset_llk(i1 = iel[[ii]], i2 = iel[[ii]], deriv = deriv),
                                       deriv = deriv))
  
  d1 <- do.call("c", lapply(der, "[[", "d1"))
  
  d2 <- NULL 
  if( deriv > 1){ 
    
    # Hessian blocks
    H <- vector(mode = "list", length = nc^2)
    for( ir in 1:nc ){
      for( ic in ir:nc ){
        Hess <- paste0(".Hess.", class(o$eff[[ic]]), "_", class(o$eff[[ir]]))
        H[[ic+(ir-1)*nc]] <- if(ir != ic){ 
          do.call(Hess, list("o1" = o$eff[[ic]], "o2" = o$eff[[ir]], 
                             "llk" = subset_llk(i1 = iel[[ic]], i2 = iel[[ir]], deriv = 2)))
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
  
  out <- list("d1" = d1, "d2" = d2)
  
  return( out )
  
}