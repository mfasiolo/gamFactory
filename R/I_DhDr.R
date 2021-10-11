
###
.DHDr <- function(o, llk, DbDr, index){
  
  neff <- length( unique(index) )
  
  # Only one effect involved: call effect-specific method
  if( neff == 1 ){ 
    llk$d2 <- llk$d2[[1]]
    DHnam <- paste0("DHessDrho.", class(o[[1]])[1])
    out <- do.call(DHnam, list(o = o[[1]], llk = llk, DbDr = DbDr))[[1]]
    
    
  } else {
    
    out <- DHessDrho.easy(o = o, llk = llk, DbDr = DbDr, index = index)
    
  }
  
  return( out )
  
}

.transEta <- function(eta, o, ii){
  if(class(o)[1] %in% c("si", "nexpsm", "mgks")){
    if(ii == 1){ eta <- eta * o$store$f1 }
  } else {
    stopifnot(class(o) == "stand")
  } 
  return( eta )
}

.getM <- function(o, ii){
  if(class(o)[1] %in% c("si", "nexpsm", "mgks")){
    M <- if(ii == 1){  o$store$g1 } else { o$store$X0 }
  } else {
    M <- o$store$X
    stopifnot(class(o) == "stand")
  } 
  return( M )
}

.getDetaDr <- function(o, DbDr, ii){
  M <- .getM(o, ii)
  if(class(o)[1] %in% c("si", "nexpsm", "mgks")){
    if(ii == 1){ 
      Deta <- M %*% DbDr[1:o$na] 
    } else {
      Deta <- M %*% DbDr[-(1:o$na)]  
    }
  } else {
    Deta <- M %*% DbDr
    stopifnot(class(o) == "stand")
  } 
  return( drop(Deta) )
}

DHessDrho.easy <- function(o, llk, DbDr, index){
  
  type <- sapply(o, function(inp) class(inp)[1])
  
  # Indexes of duplicate nested effects (same effect twice)
  duplNE <- which(type != "stand" & (index %in% index[duplicated(index)]))
  duplNENL <- which(!(type %in% c("stand", "si")) & (index %in% index[duplicated(index)]))
  notDupl <- setdiff(1:3, duplNE)
  
  # The _same_ nested appears twice
  twoNE <- length(duplNE) == 2
  twoNENL <- length(duplNENL) == 2
  
  # o[[1]] and o[[2]] are the _same_ nested effect
  needD <- twoNE && all(duplNE == c(1, 2))
  
  # o[[1]] and o[[2]] are the _same_ nested _non-linear_ effect
  needVaa3 <- twoNENL && needD 
  
  # o[[1]] and o[[3]] OR o[[2]] and o[[3]] are the _same_ nested effect
  needEG <- twoNE && !needD
  
  # o[[1]] and o[[3]] OR o[[2]] and o[[3]] are the _same_  nested _non-linear_ effect
  needMVQ <- twoNENL && needEG
  
  # Standard effects have one linear predictors, nested effects 2
  nc <- (type != "stand") + 1
  
  # Derivative wrt first-and-second and first-and-third linear predictors
  llkH12 <- llk$d2[[2]]
  llkH13 <- llk$d2[[3]]

  out <- NULL
  dH <- list()
  for(i3 in 1:nc[3]){
    jj <- 1
    V0 <- .transEta(llk$d3, o[[3]], i3) * .getDetaDr(o[[3]], DbDr, i3) 
    for(i2 in 1:nc[2]){
      V1 <- .transEta(V0, o[[2]], i2)
      M2 <- .getM(o[[2]], i2)
      for(i1 in 1:nc[1]){
        V2 <- .transEta(V1, o[[1]], i1)
        # We have two SIs we are looking at alpha entry. Hence we need additional entry (NOT YET DOCUMENTED IN TEXT)
        if( twoNE && all(c(i1, i2, i3)[duplNE] == c(1, 1)) ){
          tmp <- if(notDupl %in% c(1, 2)){ llkH12 } else{ llkH13 } 
          V2 <- V2 + .transEta(tmp * o[[duplNE[1]]]$store$f2, 
                              o[[notDupl]], c(i1, i2, i3)[notDupl]) * .getDetaDr(o[[3]], DbDr, i3)  
        }
        # Need to add E to diagonal
        if( (needEG && i3 == 2) && ( (index[3] == index[2] && i2 == 1) || (index[3] == index[1] && i1 == 1) ) ){
          V2 <- V2 + .transEta(llkH12, o[[notDupl]], c(i1, i2, i3)[notDupl]) * drop(o[[3]]$store$X1 %*% DbDr[-(1:o[[3]]$na)])
        }
        M1 <- .getM(o[[1]], i1) 
        M2VM1 <- crossprod(M2, V2 * M1)
        # Need to add G factor
        if( (needEG && i3 == 1) ){
          if(index[3] == index[2] && i2 == 2){
            D <- .transEta(llkH12, o[[notDupl]], c(i1, i2, i3)[notDupl]) * .getDetaDr(o[[3]], DbDr, i3)
            M2VM1 <- M2VM1 + crossprod(o[[2]]$store$X1, D * M1)  
          } 
          if(index[3] == index[1] && i1 == 2){
            D <- .transEta(llkH12, o[[notDupl]], c(i1, i2, i3)[notDupl]) * .getDetaDr(o[[3]], DbDr, i3)
            M2VM1 <- M2VM1 + crossprod(M2, D * o[[1]]$store$X1)
          }
        }
        if( needMVQ && i3 == 1 ) {
          if( notDupl == 1 && i2 == 1 ) {
            Q <- .contract2(X = o[[3]]$store$g2, v = DbDr[1:o[[3]]$na], ind = trind.generator(o[[3]]$na)$i2)
            D <- .transEta(.transEta(llkH13, o[[1]], i1), o[[3]], 1)
            M2VM1 <- M2VM1 + crossprod(Q, D * M1)
          } 
          if( notDupl == 2 && i1 == 1 ) {
            Q <- .contract2(X = o[[3]]$store$g2, v = DbDr[1:o[[3]]$na], ind = trind.generator(o[[3]]$na)$i2)
            D <- .transEta(.transEta(llkH12, o[[2]], i2), o[[3]], 1)
            M2VM1 <- M2VM1 + crossprod(M2, D * Q)
          }
        }
        dH[[jj]] <- if(i3 == 1) { M2VM1 } else { dH[[jj]] + M2VM1 }
        jj <- jj + 1
      }
    }
    
    # Extra X'^T %*% D %*% g1 component needed 
    if( needD ){
      D <- .transEta(llkH13, o[[3]], i3) * .getDetaDr(o[[3]], DbDr, i3)
      tmp <- crossprod(o[[1]]$store$g1, D * o[[1]]$store$X1)
      dH[[2]] <- dH[[2]] + tmp
      dH[[3]] <- dH[[3]] + t(tmp)
    }
    if( needVaa3 ){
      D <- .transEta(.transEta(llkH13, o[[1]], 1), o[[3]], i3) * .getDetaDr(o[[3]], DbDr, i3)
      tmp <- .vec_to_sym_mat(colSums(D * o[[1]]$store$g2), o[[1]]$na)
      dH[[1]] <- dH[[1]] + tmp
    }
  }
  
  # Build Hessian # INEFFICIENT
  for(ir in 1:nc[2]){
    out <- rbind(out, do.call("cbind", dH[(1+(ir-1)*nc[1]):(ir*nc[1])]))
  }
  
  return( out )
  
}

