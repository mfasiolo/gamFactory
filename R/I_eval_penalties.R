########
# Evaluate additional penalties and their derivatives
#
.eval_penalties <- function(eff, info, d1b, deriv, outer){
  effType <- sapply(info$type, paste0, collapse = '')
  nstand <- which(effType != "stand")
  pen <- vector(mode = "list", length = length(nstand))
  kk <- 1
  for(ii in nstand){
    iec <- info$iec[[ii]]
    extra <- info$extra[[ii]]
    aii <- iec[1:eff[[ii]]$na]
    pen[[kk]] <- pen_var(o = eff[[ii]], v = 1, deriv = deriv)
    if(outer){
      pen[[kk]]$outer <- pen_var_outer(o = eff[[ii]], v = 1, DaDr = d1b[aii, , drop = FALSE])
    }
    pen[[kk]]$iec <- aii 
    kk <- kk + 1
  }
  return( pen )
}
