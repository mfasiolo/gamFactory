#'
#' Derivative of Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.linearPredictor
#' @rdname DHessDrho.linearPredictor
#' @export DHessDrho.linearPredictor
#'
DHessDrho.linearPredictor <- function(o, llk, DbDr){
  
  der <- o$xtra$der
  KK <- o$xtra$KK
  
  # Number of components and smoothing parameters
  nc <- length( der )
  m <- ncol( DbDr )
  
  # Derivatives of Hessian blocks wrt smoothing param 
  # INEFFICIENT should calculate only upper triangle
  out <- lapply(1:m, function(nouse) NULL)
  for( ism in 1:m ){ # Loop over smoothing parameters
    jj <- 1
    dH <- list()
    for( ir in 1:nc ){
      for( ic in 1:nc ){
        for( i3 in 1:nc ){
          tmp <- .DHDr(o = der[ c(ic, ir, i3) ],
                       llk = llk, 
                       DbDr = DbDr[KK[[i3]], ism, drop = FALSE], 
                       index = c(ic, ir, i3))
          dH[[jj]] <- if( i3 == 1 ){ tmp } else { dH[[jj]] + tmp }
        }
        jj <- jj + 1
      }
    }
    
    # Build Hessian # INEFFICIENT
    for(ir in 1:nc){
      out[[ism]] <- rbind(out[[ism]], do.call("cbind", dH[(1+(ir-1)*nc):(ir*nc)]))
    }
    
  }
  
  return( out )
  
}