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
    if( "nested" %in%  info$type[[ii]] ){
      
      ism <- info$extra[[ii]]$ism
      sii <- o$smooth[[ism]]
      sii$xt$si$alpha <- coef(o)[ info$iec[[ii]][1:length(sii$xt$si$alpha)] ]
      
      # Inner smooth must be centered using original data
      if(is.null(sii$xt$si$xm)){
        sii$xt$si$xm <- mean(attr(Predict.matrix.nested(sii, data = o$model), "inner_linpred_unscaled"))
      }
      
      o$smooth[[ism]] <- sii
      
   }
    
  }
  
  o$linear.predictors <- predict(o, type = "link")
  o$fitted.values <- predict(o, type = "response")
  
  return( o )
  
}