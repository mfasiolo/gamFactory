#' 
#' Post-process a GAM containing non-standard effects
#' @name postproc_gam_nl
#' @rdname postproc_gam_nl
#' @export postproc_gam_nl
#'
postproc_gam_nl <- function(o, info){
  
  sms <- o$smooth
  ne <- length( info$type )
  for(ii in 1:ne){
    
    # Coefficients of nested smooth must be saved in smooth object
    if( ("si.smooth" %in%  info$type[ii]) || 
        ("nexpsm.smooth" %in%  info$type[ii]) ||
        ("mgks.smooth" %in%  info$type[ii]) ){
      ism <- info$extra[[ii]]$ism
      sms[[ism]]$xt$si$alpha <- coef(o)[ info$iec[[ii]][1:length(sms[[ism]]$xt$si$alpha)] ]
    }
    
  }
  
  o$smooth <- sms
  
  o$fitted.values <- predict(o, type = "response")
  
  return( o )
  
}