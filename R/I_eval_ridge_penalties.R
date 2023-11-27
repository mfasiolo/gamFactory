########
# Evaluate additional ridge penalties and their derivatives
#
.eval_ridge_penalties <- function(eff, info, deriv){
  effType <- sapply(info$type, paste0, collapse = '')
  nstand <- which(effType != "stand")
  pen <- vector(mode = "list", length = length(nstand))
  kk <- 1
  for(ii in nstand){
    iec <- info$iec[[ii]]
    extra <- info$extra[[ii]]
    aii <- iec[1:eff[[ii]]$na]
    cl <- class(eff[[ii]])
    if("nested" %in% cl){
      if("si" %in% cl){
        ipc <- 1:eff[[ii]]$na
      } else{
        ipc <- 2:eff[[ii]]$na
      }
        pen[[kk]] <- pen_ridge_var(o = eff[[ii]], extra = extra, ipc = ipc, deriv = deriv)
    }
    
    pen[[kk]]$iec <- aii 
    kk <- kk + 1
  }
  return( pen )
}

pen_ridge_var <- function(o, extra, ipc, deriv){
  
  a <- o$param[ipc]
  a_init <- extra$si$alpha[ipc]
  
  l0 <- .5 * sum((a-a_init)^2)
  
  l1 <- l2 <- NULL
  if(deriv){
    l1 <- rep(0, o$na)
    l1[ipc] <- a - a_init
    
    if(deriv > 1){
      l2 <- matrix(0, o$na, o$na)
      l2[ipc,ipc] <- 1
    }
  }
  return(list("d0" = l0, "d1" = l1, "d2" = l2))
}
