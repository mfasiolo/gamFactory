
.wrap_nested_basis <- function(b, P, slope){
  
  force(b); force(P); force(slope)
  
  evalX <- function(x, deriv){
    withCallingHandlers({
      o <- b$evalX(x = x, deriv = deriv)
      if(slope){ o$X0 <- cbind(x, o$X0) }
      o$X0 <- o$X0 %*% P
      if(deriv > 0){
        if(slope){ o$X1 <- cbind(1, o$X1) }
        o$X1 <- o$X1 %*% P
        if(deriv > 1){
          if(slope){ o$X2 <- cbind(0, o$X2) }
          o$X2 <- o$X2 %*% P
        }
        if(deriv > 2){
          if(slope){ o$X3 <- cbind(0, o$X3) }
          o$X3 <- o$X3 %*% P
        }
      }
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