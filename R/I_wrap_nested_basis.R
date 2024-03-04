
.wrap_nested_basis <- function(b, P, Xth, add_slope){
  
  force(b); force(P); force(Xth)
  
  evalX <- function(x, deriv){
    withCallingHandlers({
      o <- b$evalX(x = x, deriv = deriv) # Get raw basis & derivatives
      if(add_slope){ # Add a slope in the last column: 1st derivative is 1, rest 0.
        o$X0 <- cbind(o$X0, x)
        if(deriv){
          o$X1 <- cbind(o$X1, 1)
          if(deriv >= 2){
            for(ii in 2:deriv){
              o[[paste0("X", ii)]] <- cbind(o[[paste0("X", ii)]], 0)
            }
          }
        }
      }
      o <- lapply(o, function(X) X %*% P)       # Reparametrise
      o <- linextr(x = x, b = o, th = b$krange, # Linearly extrapolate: NOTE method = "smooth" won't work if add_slope == TRUE
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