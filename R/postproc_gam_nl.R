#' 
#' Post-process a GAM containing non-standard effects
#' @name postproc_gam_nl
#' @rdname postproc_gam_nl
#' @export postproc_gam_nl
#'
postproc_gam_nl <- function(o, info){
  
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
      
      if("si_nexpsm" %in%  info$type[[ii]]){
        n_si <- sii$xt$si$n_si
        n_nexp <- sii$xt$si$n_nexp
        sii$xt$si$alpha_nexp <- as.vector(sii$xt$si$alpha[1:n_nexp])
        sii$xt$si$alpha_si <- as.vector(sii$xt$si$alpha[(n_nexp+1):(n_si+n_nexp)])
        
        if (isTRUE(sii$xt$si$positive_si)) {
          sii$xt$si$alpha_si_inner <- sii$xt$si$alpha_si        # optimized alpha_si
          sii$xt$si$alpha_si_true  <- exp(sii$xt$si$alpha_si)   # positive alpha_si
        }
        
        sii$xt$si$xm <- c(sii$xt$si$xm, 
                          nexp = (mean(attr(Predict.matrix.nested(sii, data = o$model), 
                                    "inner_linpred_unscaled"))))
      }
      
      jacobian <- get_jacobian.nested(sii, data = o$model, param = coef(o)[info$iec[[ii]]])
      sii$xt$jacobian <- jacobian$JJ
      sii$xt$xa <- jacobian$xa  #smoothed s_t, useful for f(s_t) vs s_t plot
      
      o$smooth[[ism]] <- sii
      
    }
    
  }

  o$linear.predictors <- predict.gamnl(o, type = "link")
  o$fitted.values <- predict.gamnl(o, type = "response")
  
  return( o )
  
}
