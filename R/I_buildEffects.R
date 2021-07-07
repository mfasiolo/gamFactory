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
      
      eff[[ii]] <- eff_stand( X[ , iec, drop = FALSE] )  
      
    } else {
      
      extra <- effInfo$extra[[ii]]
      Xi <- extra$si$X
      vr <- extra$si$vr
      
      if( effType[ii] == "si.smooth" ){
        
        aii <- iec[ 1:ncol(Xi) ]
        eff[[ii]] <- eff_si(Xi = Xi, splineDes = extra$splineDes) 
        pen[[kk]] <- pen_varSI(a = coef[aii], x = Xi, v = vr, deriv = deriv)
        if(outer){
          pen[[kk]]$outer <- pen_varSI_outer(a = coef[aii], x = Xi, DaDr = d1b[aii, , drop = FALSE])
        }
        pen[[kk]]$iec <- aii
        
        kk <- kk + 1
        
      } else {
        
        if( effType[ii] == "nexpsm.smooth" ){
          
          eff[[ii]] <- eff_nexpsm(y = extra$si$x, Xi = Xi, splineDes = extra$splineDes) 
   
        } else {
          
          stop("Don't know this effect type")  
          
        }
        
      }
    }
  }
  
  return( list("eff" = eff, "pen" = pen) )
  
}