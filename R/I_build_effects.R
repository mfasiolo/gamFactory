############################
############################
# Returns list of effects 
#
.build_effects <- function(X, info, outer){
  effType <- sapply(info$type, paste0, collapse = '.')
  ne <- length( effType )
  eff <- list( )
  for(ii in 1:ne){
    iec <- info$iec[[ii]]
    if( effType[ii] == "stand" ){
      eff[[ii]] <- eff_stand( X[ , iec, drop = FALSE] )  
    } else {
      extra <- info$extra[[ii]]
      Xi <- extra$si$X
      if( effType[ii] == "si.nested" ){
        eff[[ii]] <- eff_si(Xi = Xi, basis = extra$basis) 
      } else {
        if( effType[ii] == "nexpsm.nested" ){
          eff[[ii]] <- eff_nexpsm(y = extra$si$x, Xi = Xi, basis = extra$basis, times = extra$si$times) 
        } else {
          if( effType[ii] == "mgks.nested" ){
            eff[[ii]] <- eff_mgks(y = extra$si$x, dist = extra$si$dist, basis = extra$basis)
          } else {
          stop("Don't know this effect type")  
          }
        }
      }
    }
  }
  return( eff )
}

