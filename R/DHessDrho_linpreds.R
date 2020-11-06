#'
#' Derivative of Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.linpreds
#' @rdname DHessDrho.linpreds
#' @export 
#'
DHessDrho.linpreds <- function(o, llk, DbDr){
  
  eff <- o$eff
  nlp <- o$nlp
  iec <- o$iec
  iel <- o$iel
  nc  <- o$nc
  
  ord <- order(do.call("c", iec))
  
  # Number of smoothing parameters
  m <- ncol( DbDr )
  
  subset_llk <- .create_subset_llk(llk, nlp) 
  
  # Derivatives of Hessian blocks wrt smoothing param 
  out <- lapply(1:m, function(nouse) NULL)
  for( ism in 1:m ){ # Loop over smoothing parameters
    jj <- 1
    dH <- vector(mode = "list", length = nc^2)
    for( ir in 1:nc ){
      for( ic in ir:nc ){
        for( iacc in 1:nc ){
          tmp <- .DHDr(o = eff[ c(ic, ir, iacc) ],
                       llk = subset_llk(i1 = iel[[ic]], i2 = iel[[ir]], i3 = iel[[iacc]], deriv = 3), 
                       DbDr = DbDr[iec[[iacc]], ism, drop = FALSE], 
                       index = c(ic, ir, iacc))
          dH[[ic+(ir-1)*nc]] <- if( iacc == 1 ){ tmp } else { dH[[ic+(ir-1)*nc]] + tmp }
        }
        if( ir != ic ){
          dH[[ir+(ic-1)*nc]] <- t( dH[[ic+(ir-1)*nc]] )
        }
        jj <- jj + 1
      }
    }
    
    # Build Hessian # INEFFICIENT
    for(ir in 1:nc){
      out[[ism]] <- rbind(out[[ism]], do.call("cbind", dH[(1+(ir-1)*nc):(ir*nc)]))
    }
    
    out[[ism]] <- out[[ism]][ord, ord]
    
  }

  return( out )
  
}