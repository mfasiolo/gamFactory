
.wrap_nested_basis <- function(b, P, Xth){
  
  force(b); force(P); force(Xth)
  
  evalX <- function(x, deriv){
    withCallingHandlers({
      o <- b$evalX(x = x, deriv = deriv)        # Get raw basis & derivatives
      o <- lapply(o, function(X) X %*% P)       # Reparametrise
      o <- linextr(x = x, b = o, th = b$krange, # Linearly extrapolate
                   Xbo = Xth$X0%*%P, Xbo1 = Xth$X1%*%P, method = "simple")
      o
    }, warning = function(w) {
      if (length(grep("there is \\*no\\* information about some basis coefficients", conditionMessage(w)))){
        invokeRestart("muffleWarning")
      }
    })
  }
  
  out <- list("evalX" = evalX)

  return(out)
}