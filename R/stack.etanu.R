# This function converts le lee leee to ln lnn lnnn when using stack effects
# K is the number of parameters nu.
stack.etanu <- function(le = NULL, lee = NULL, leee = NULL, 
                        en = NULL, enn = NULL, ennn = NULL,
                        deriv = 1) {
  
  out <- list()
  
  if (1 %in% deriv) {
    out$ln <- le * en
  }
  
  if (2 %in% deriv) {
    i2 <- trind.generator(K - 1)$i2
    
    lnn <- list()
    coun <- 1
    for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
      lnn[[coun]] <- lee * en[, rr] * en[, ss] + le * enn[, i2[rr, ss]]
      coun <- coun + 1
    }
    out$lnn <- do.call(cbind, lnn)
  }
  
  if (3 %in% deriv) {
    i2 <- trind.generator(K - 1)$i2
    i3 <- trind.generator(K - 1)$i3
    
    lnnn <- list()
    coun <- 1
    for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) for (tt in ss:(K - 1)) {
      lnnn[[coun]] <- leee * en[, rr] * en[, ss] * en[, tt] +
        lee * (enn[, i2[rr, tt]] * en[, ss] + 
                 enn[, i2[ss, tt]] * en[, rr] + 
                 enn[, i2[rr, ss]] * en[, tt]) +
        le * ennn[, i3[rr, ss, tt]]
      coun <- coun + 1
    }
    out$lnnn <- do.call(cbind, lnnn)
  }
  
  return(out)
}

