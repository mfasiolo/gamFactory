############################
############################
# Returns list of effects and penalties
#
#
.buildEffects <- function(X, coef, effInfo, d1b, deriv, outer){
  
  effType <- effInfo$type
  ne <- length( effType )
  
  eff <- pen <- list()
  
  kk <- 1
  # Loop through all the effects
  for(ii in 1:ne){
    
    iec <- effInfo$iec[[ii]]
    
    if( effType[ii] == "standard" ){
      
      eff[[ii]] <- buildStandardEffect( X[ , iec, drop = FALSE] )  
      
    } else {
      
      extra <- effInfo$extra[[ii]]
      
      if( effType[ii] == "singleIndex" ){
        
        aii <- ncol(extra$Xi)
        eff[[ii]] <- buildSingleIndexEffect(Xi = extra$Xi, 
                                            splineDes = extra$splineDes) 
        pen[[kk]] <- pen_varSI(a = coef[aii], 
                               x = extra$Xi, 
                               v = extra$vr, 
                               deriv = deriv)
        if(outer){
          pen[[kk]]$outer <- pen_varSI_outer(a = coef[aii], 
                                             x = extra$Xi, 
                                             DaDr = d1b[aii, , drop = FALSE])
        }
        pen[[kk]]$iec <- aii
        
      } else {
        
        stop("Don't know this effect type")
        
      }
    }
  }
  
  return( list("eff" = eff, "pen" = pen) )
  
}