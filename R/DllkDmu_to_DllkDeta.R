######
# Chain rule: transform derivatives w.r.t. mu to derivatives w.r.t. eta
######
DllkDMu_to_DllkDeta <- function(DllkDMu, etas, mus, family, wt, deriv){
  
  linfo <- family$linfo
  out <- NULL
  np <- length( mus )
  
  if (deriv >= 1) {
    
    l2 <- l3 <- l4 <- g2 <- g3 <- g4 <- i2 <- i3 <- i4 <- 0 ## defaults
    
    l1 <- wt  * do.call("cbind", DllkDMu$d1)
    ig1 <- do.call("cbind", lapply(1:np, function(.kk) linfo[[.kk]]$mu.eta( etas[[.kk]])))
    
    if( deriv >= 2 ){
      
      l2 <- wt  * do.call("cbind", DllkDMu$d2)
      g2  <- do.call("cbind", lapply(1:np, function(.kk) linfo[[.kk]]$d2link( mus[[.kk]])))
      i2 <- family$tri$i2 
      
      if( deriv >= 3 ){
        
        l3 <- wt  * do.call("cbind", DllkDMu$d3) 
        g3  <- do.call("cbind", lapply(1:np, function(.kk) linfo[[.kk]]$d3link( mus[[.kk]])))
        i3 <- family$tri$i3
        
        if( deriv >= 4 ){
          
          l4 <- wt  * do.call("cbind", DllkDMu$d4) 
          g4 <- do.call("cbind", lapply(1:np, function(.kk) linfo$linfo[[.kk]]$d4link( mus[[.kk]])))
          i4 <- family$tri$i4
          
        }
        
      }
      
    }
    
    out <- gamlss.etamu(l1, l2, l3, l4, ig1, g2, g3, g4, i2, i3, i4, 
                        deriv = switch(as.character(deriv), "0" = 0, "1" = 0, "2" = 0, "3" = 1, "4" = 3))
    
    if( is.matrix(out$l1) ) out$d1 <- lapply(seq_len(ncol(out$l1)), function(ii) out$l1[ , ii])
    if( is.matrix(out$l2) ) out$d2 <- lapply(seq_len(ncol(out$l2)), function(ii) out$l2[ , ii])
    if( is.matrix(out$l3) ) out$d3 <- lapply(seq_len(ncol(out$l3)), function(ii) out$l3[ , ii])
    if( is.matrix(out$l4) ) out$d4 <- lapply(seq_len(ncol(out$l4)), function(ii) out$l4[ , ii])
    
  }
  
  
  return( out )
  
}