
###
.DHDr <- function(o, llk, DbDr, index){
  
  neff <- length( unique(index) )
  
  # Only one effect involved: call effect-specific method
  if( neff == 1 ){ 
    llk$d2 <- llk$d2[[1]]
    DHnam <- paste0("DHessDrho.", class(o[[1]]))
    out <- do.call(DHnam, list(o = o[[1]], llk = llk, DbDr = DbDr))[[1]]
    
    
  } else {
    
    out <- DHessDrho.easy(o = o, llk = llk, DbDr = DbDr, index = index)
    
  }
  
  return( out )
  
}

transEta <- function(eta, o, ii){
  if(class(o) == "singleIndex"){
    if(ii == 1){ eta <- eta * o$store$f1 }
  } else {
    stopifnot(class(o) == "standard")
  } 
  return( eta )
}

getM <- function(o, ii){
  if(class(o) == "singleIndex"){
    M <- if(ii == 1){  o$store$Xi } else { o$store$X0 }
  } else {
    M <- o$store$X
    stopifnot(class(o) == "standard")
  } 
  return( M )
}

getDetaDr <- function(o, DbDr, ii){
  M <- getM(o, ii)
  if(class(o) == "singleIndex"){
    if(ii == 1){ 
      Deta <- M %*% DbDr[1:o$na] 
    } else {
      Deta <- M %*% DbDr[-(1:o$na)]  
    }
  } else {
    Deta <- M %*% DbDr
    stopifnot(class(o) == "standard")
  } 
  return( drop(Deta) )
}

DHessDrho.easy <- function(o, llk, DbDr, index){
  
  type <- sapply(o, class)
  
  # Indexes of duplicate single index effect (same effect twice)
  duplSI <- which(type == "singleIndex" & (index %in% index[duplicated(index)]))
  notDupl <- setdiff(1:3, duplSI)
  
  # The _same_ single index appears twice
  twoSI <- length(duplSI) == 2
  
  # o[[1]] and o[[2]] are the _same_ single index effect
  needD <- twoSI && all(duplSI == c(1, 2))
  
  # o[[1]] and o[[3]] OR o[[2]] and o[[3]] are the _same_ single index effect
  needEG <- twoSI && !needD
  
  # Standard effects have one linear predictors, nested effects 2
  nc <- (type != "standard") + 1
  
  # Derivative wrt first-and-second and first-and-third linear predictors
  llkH12 <- llk$d2[[2]]
  llkH13 <- llk$d2[[3]]

  out <- NULL
  dH <- list()
  for(i3 in 1:nc[3]){
    jj <- 1
    V0 <- transEta(llk$d3, o[[3]], i3) * getDetaDr(o[[3]], DbDr, i3) 
    for(i2 in 1:nc[2]){
      V1 <- transEta(V0, o[[2]], i2)
      M2 <- getM(o[[2]], i2)
      for(i1 in 1:nc[1]){
        V2 <- transEta(V1, o[[1]], i1)
        # We have two SIs we are looking at alpha entry. Hence we need additional entry (NOT YET DOCUMENTED IN TEXT)
        if( twoSI && all(c(i1, i2, i3)[duplSI] == c(1, 1)) ){
          tmp <- if(notDupl %in% c(1, 2)){ llkH12 } else{ llkH13 } 
          V2 <- V2 + transEta(tmp * o[[duplSI[1]]]$store$f2, 
                              o[[notDupl]], c(i1, i2, i3)[notDupl]) * getDetaDr(o[[3]], DbDr, i3)  
        }
        # Need to add E to diagonal
        if( (needEG && i3 == 2) && ( (index[3] == index[2] && i2 == 1) || (index[3] == index[1] && i1 == 1) ) ){
          V2 <- V2 + transEta(llkH12, o[[notDupl]], c(i1, i2, i3)[notDupl]) * drop(o[[3]]$store$X1 %*% DbDr[-(1:o[[3]]$na)])
        }
        M1 <- getM(o[[1]], i1) 
        M2VM1 <- crossprod(M2, V2 * M1)
        # Need to add G factor
        if( (needEG && i3 == 1) ){
          if(index[3] == index[2] && i2 == 2){
            D <- transEta(llkH12, o[[notDupl]], c(i1, i2, i3)[notDupl]) * getDetaDr(o[[3]], DbDr, i3)
            M2VM1 <- M2VM1 + crossprod(o[[2]]$store$X1, D * M1)  
          } 
          if(index[3] == index[1] && i1 == 2){
            D <- transEta(llkH12, o[[notDupl]], c(i1, i2, i3)[notDupl]) * getDetaDr(o[[3]], DbDr, i3)
            M2VM1 <- M2VM1 + crossprod(M2, D * o[[1]]$store$X1)
          }
        }
        dH[[jj]] <- if(i3 == 1) { M2VM1 } else { dH[[jj]] + M2VM1 }
        jj <- jj + 1
      }
    }
    
    # Extra X'^T %*% D %*% Xi component needed 
    if( needD ){
      D <- transEta(llkH13, o[[3]], i3) * getDetaDr(o[[3]], DbDr, i3)
      tmp <- crossprod(o[[1]]$store$Xi, D * o[[1]]$store$X1)
      dH[[2]] <- dH[[2]] + tmp
      dH[[3]] <- dH[[3]] + t(tmp)
    }
  }
  
  # Build Hessian # INEFFICIENT
  for(ir in 1:nc[2]){
    out <- rbind(out, do.call("cbind", dH[(1+(ir-1)*nc[1]):(ir*nc[1])]))
  }
  
  return( out )
  
}

