###################################
# Extract the part of the log-likelihood relevant for (mixed) derivatives wrt linear predictors i1, i2, i3
.create_subset_llk <- function(llk, p){
  
  indexMat <- trind.generator(p)
  
  subset_llk <- function(i1, i2 = NULL, i3 = NULL, deriv = 1){
    if( !is.list(llk$d1) ){ return(llk) }
    if( deriv >= 1 ){
      llk$d1 <- llk$d1[[ i1 ]]
      if( deriv >= 2 ){
        llk$d2 <- llk$d2[[ indexMat$i2[i1, i2] ]]
        if( deriv >= 3){
          llk$d3 <- llk$d3[[ indexMat$i3[i1, i2, i3] ]]
        }
      }
    }
    return(llk)
  }
  
  return( subset_llk ) 
  
}

