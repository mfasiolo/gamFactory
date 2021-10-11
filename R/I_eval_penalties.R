########
# Evaluate additional penalties and their derivatives
#
.eval_penalties <- function(eff, info, d1b, deriv, outer){
  effType <- info$type
  nstand <- which(info$type != "stand")
  pen <- vector(mode = "list", length = length(nstand))
  kk <- 1
  for(ii in nstand){
    iec <- info$iec[[ii]]
    extra <- info$extra[[ii]]
    vr <- extra$si$vr
    aii <- iec[1:eff[[ii]]$na]
    pen[[kk]] <- pen_var(o = eff[[ii]], v = vr, deriv = deriv)
    if(outer){
      pen[[kk]]$outer <- pen_var_outer(o = eff[[ii]], v = vr, DaDr = d1b[aii, , drop = FALSE])
    }
    pen[[kk]]$iec <- aii 
    kk <- kk + 1
  }
  return( pen )
}
