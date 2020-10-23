.postProcNonLinGAM <- function(o, info){
  
  sms <- o$smooth
  ne <- length( info$type )
  for(ii in 1:ne){
    
    # Coefficients of single index must be saved in smooth object
    if( "si.smooth" %in%  info$type[ii] ){
      ism <- info$extra[[ii]]$ism
      sms[[ism]]$xt$si$alpha <- coef(o)[ info$iec[[ii]][1:ncol(sms[[ism]]$xt$si$X)] ]
    }
    
  }
  
  o$smooth <- sms
  
  o$fitted.values <- predict(o, type = "response")
  
  return( o )
  
}