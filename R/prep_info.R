#' 
#' Get information to build family for GAM with non-standard effect 
#' @name prep_info
#' @rdname prep_info
#' @export prep_info
#' 
prep_info <- function(o){
  
  lpi <- attr(o$X, "lpi")
  
  # Coefficients to effects correspondence
  if( !is.list(lpi) ){ lpi <- list(lpi) }
  nlp <- length( lpi )
  sms <- o$smooth
  cls <- lapply(sms, class)
  sp <- which(sapply(sms, function(.x){ return(!is.null(.x$xt$si)) }))
  nsp <- length(sp)
  spCo <- lapply(sp, function(.kk) sms[[.kk]]$first.para:sms[[.kk]]$last.para) 
  
  out <- list("iec" = lpi, 
              "iel" = c(1:nlp, rep(NA, nsp)), 
              "type" = as.list(c(rep("stand", nlp), rep(NA, nsp))), 
              "extra" = rep( list(list()), nlp + nsp ))
  
  if( nsp ){
    for( ii in 1:nsp ){
      # Index of lin pred to which ii-th special effect belongs
      kk <- which( sapply(lpi, function(.ff) any(.ff %in% spCo[[ii]]) ) )
      out$iec[[kk]] <- out$iec[[kk]][ !(out$iec[[kk]] %in% spCo[[ii]]) ]  
      out$iec[[nlp+ii]] <- spCo[[ii]]
      out$iel[nlp+ii] <- kk
      out$type[[nlp+ii]] <- cls[[ sp[ii] ]][1:2]
      out$extra[[nlp+ii]] <- sms[[ sp[ii] ]]$xt
      out$extra[[nlp+ii]]$ism <- sp[[ ii ]] 
    }
  }
  return(out)
}