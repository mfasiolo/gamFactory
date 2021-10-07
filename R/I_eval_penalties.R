########
# Evaluate additional penalties and their derivatives
#
.eval_penalties <- function(eff, effInfo, d1b, deriv, outer){
  effType <- effInfo$type
  nstand <- which(effInfo$type != "standard")
  pen <- vector(mode = "list", length = length(nstand))
  kk <- 1
  for(ii in nstand){
    iec <- effInfo$iec[[ii]]
    extra <- effInfo$extra[[ii]]
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
