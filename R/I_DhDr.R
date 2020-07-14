
###
.DHDr <- function(o, llk, DbDr, index){
  
  neff <- length( unique(index) )
  
  # Only one effect involved: call effect-specific method
  if( neff == 1 ){ 
    
    DHnam <- paste0("DHessDrho.", o[[1]]$type)
    out <- do.call(DHnam, list(o = o[[1]], llk = llk, DbDr = DbDr))[[1]]
    
    
  } else {
    
    out <- DHessDrho.easy(o = o, llk = llk, DbDr = DbDr, index = index)
    
  }
  
  return( out )
  
}

transEta <- function(eta, o, ii){
  if(o$type == "singleIndex"){
    if(ii == 1){ eta <- eta * o$xtra$f1 }
  } else {
    stopifnot(o$type == "standard")
  } 
  return( eta )
}

getM <- function(o, ii){
  if(o$type == "singleIndex"){
    M <- if(ii == 1){  o$xtra$Xi } else { o$xtra$X }
  } else {
    M <- o$xtra$X
    stopifnot(o$type == "standard")
  } 
  return( M )
}

getDetaDr <- function(o, DbDr, ii){
  M <- getM(o, ii)
  if(o$type == "singleIndex"){
    if(ii == 1){ 
      Deta <- M %*% DbDr[1:o$xtra$na] 
    } else {
      Deta <- M %*% DbDr[-(1:o$xtra$na)]  
    }
  } else {
    Deta <- M %*% DbDr
    stopifnot(o$type == "standard")
  } 
  return( drop(Deta) )
}

DHessDrho.easy <- function(o, llk, DbDr, index){
  
  type <- sapply(o, "[[", "type")
  
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
        if( twoSI && all(c(i1, i2, i3)[duplSI] == c(1, 1)) ){
          V2 <- V2 + transEta(llk$d2 * o[[duplSI[1]]]$xtra$f2, 
                              o[[notDupl]], c(i1, i2, i3)[notDupl]) * getDetaDr(o[[3]], DbDr, i3)  
        }
        if( (needEG && i3 == 2) && ( (index[3] == index[2] && i2 == 1) || (index[3] == index[1] && i1 == 1) ) ){
          V2 <- V2 + transEta(llk$d2, o[[notDupl]], c(i1, i2, i3)[notDupl]) * drop(o[[3]]$xtra$X1 %*% DbDr[-(1:o[[3]]$xtra$na)])
        }
        M1 <- getM(o[[1]], i1) 
        M2VM1 <- crossprod(M2, V2 * M1)
        if( (needEG && i3 == 1) ){
          if(index[3] == index[2] && i2 == 2){
            D <- transEta(llk$d2, o[[notDupl]], c(i1, i2, i3)[notDupl]) * getDetaDr(o[[3]], DbDr, i3)
            M2VM1 <- M2VM1 + crossprod(o[[2]]$xtra$X1, D * M1)  
          } 
          if(index[3] == index[1] && i1 == 2){
            D <- transEta(llk$d2, o[[notDupl]], c(i1, i2, i3)[notDupl]) * getDetaDr(o[[3]], DbDr, i3)
            M2VM1 <- M2VM1 + crossprod(M2, D * o[[1]]$xtra$X1)
          }
        }
        dH[[jj]] <- if(i3 == 1) { M2VM1 } else { dH[[jj]] + M2VM1 }
        jj <- jj + 1
      }
    }
    
    # Extra X'^T %*% D %*% Xi component needed 
    if( needD ){
      D <- transEta(llk$d2, o[[3]], i3) * getDetaDr(o[[3]], DbDr, i3)
      tmp <- crossprod(o[[1]]$xtra$Xi, D * o[[1]]$xtra$X1)
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




###

##### Exceptions
# singleIndex1 singleIndex1 ...
# X V tilde(X) + X' D tilde(X) 
# alpha beta  gamma   DONE
# alpha beta  alpha*  DONE
# alpha beta  beta*   DONE

# (index[1] != index[2]) && (sum(type == "singleIndex") >= 1)



# Mixed index of same effect finishing with beta -> E
# if( (needEG && i3 == 2) && ( (index[3] == index[2] && i3 != i2) || (index[3] == index[1] && i3 != i1) ) ){
#   V2 <- V2 + transEta(llk$d2, o[[2]], i2) * getDetaDr(o[[3]], DbDr, i3)
# }



# singleIndex1 ... singleIndex1or2
# M (V + E) M
# alpha alpha* beta
# alpha alpha* beta*
# alpha beta*  beta
# alpha gamma  beta

# singleIndex1 ... singleIndex1or2
# M V M + G 
# alpha beta*  alpha*
# beta  beta*  alpha
# beta  beta*  alpha*
# beta  gamma  alpha





# DHessDrho.messy <- function(o, llk, DbDr, index){
#   
#   type <- sapply(o, "[[", "type")
#   sii <- which(type != "standard")
#   nsi <- length( sii )
#   # Standard effects have one linear predictors, nested effects 2
#   nc <- (type != "standard") + 1
#   
#   out <- NULL
#   dH <- list()
#   for(i3 in 1:nc[3]){
#     jj <- 1
#     V0 <- transEta(llk$d3, o[[3]], i3) * getDetaDr(o[[3]], DbDr, i3) 
#     for(i2 in 1:nc[2]){
#       V1 <- transEta(V0, o[[2]], i2)
#       M2 <- getM(o[[2]], i2)
#       for(i1 in 1:nc[1]){
#         V2 <- transEta(V1, o[[1]], i1)
#         M1 <- getM(o[[1]], i1) 
#         M2VM1 <- crossprod(M2, V2 * M1)
#         dH[[jj]] <- if(i3 == 1) { M2VM1 } else { dH[[jj]] + M2VM1 }
#         jj <- jj + 1
#       }
#     }
#   }
#   
#   # Build Hessian # INEFFICIENT
#   for(ir in 1:nc[2]){
#     out <- rbind(out, do.call("cbind", dH[(1+(ir-1)*nc[1]):(ir*nc[1])])) * 0
#   }
#   
#   return( out )
#   
# }


