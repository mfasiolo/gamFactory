############################
############################
# Returns list of effects 
#
.build_effects <- function(X, info, outer){
  effType <- info$type
  ne <- length( effType )
  eff <- list( )
  for(ii in 1:ne){
    iec <- info$iec[[ii]]
    if( effType[ii] == "stand" ){
      eff[[ii]] <- eff_stand( X[ , iec, drop = FALSE] )  
    } else {
      extra <- info$extra[[ii]]
      Xi <- extra$si$X
      if( effType[ii] == "si.smooth" ){
        eff[[ii]] <- eff_si(Xi = Xi, basis = extra$basis) 
      } else {
        if( effType[ii] == "nexpsm.smooth" ){
          eff[[ii]] <- eff_nexpsm(y = extra$si$x, Xi = Xi, basis = extra$basis) 
        } else {
          if( effType[ii] == "mgks.smooth" ){
            eff[[ii]] <- eff_mgks(y = extra$si$x, X = Xi, Xi = extra$si$X0, 
                                  basis = extra$basis)
          } else {
          stop("Don't know this effect type")  
          }
        }
      }
    }
  }
  return( eff )
}

