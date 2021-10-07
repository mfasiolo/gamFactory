############################
############################
# Returns list of effects 
#
.build_effects <- function(X, effInfo, outer){
  effType <- effInfo$type
  ne <- length( effType )
  eff <- list( )
  for(ii in 1:ne){
    iec <- effInfo$iec[[ii]]
    if( effType[ii] == "standard" ){
      eff[[ii]] <- eff_stand( X[ , iec, drop = FALSE] )  
    } else {
      extra <- effInfo$extra[[ii]]
      Xi <- extra$si$X
      if( effType[ii] == "si.smooth" ){
        eff[[ii]] <- eff_si(Xi = Xi, splineDes = extra$splineDes) 
      } else {
        if( effType[ii] == "nexpsm.smooth" ){
          eff[[ii]] <- eff_nexpsm(y = extra$si$x, Xi = Xi, splineDes = extra$splineDes) 
        } else {
          if( effType[ii] == "mgks.smooth" ){
            eff[[ii]] <- eff_mgks(y = extra$si$x, X = Xi, Xi = extra$si$X0, 
                                  splineDes = extra$splineDes)
          } else {
          stop("Don't know this effect type")  
          }
        }
      }
    }
  }
  return( eff )
}

