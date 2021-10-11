#'
#' Comparing analytic with numerical finite differences
#' 
#' @description Function for getting higher derivatives by finite differences.
#' @param obj the function to be differentiated.
#' @param param vector of parameters at which we should calculate the derivatives.
#' @param ord the order to derivatives to be checked.
#' @name check_deriv
#' @rdname check_deriv
#' @export check_deriv
#' @importFrom numDeriv jacobian
#'
check_deriv <- function(obj, param, ord = 1){
  
  # Helper function needed to fix some of the parameters, in order evaluate 
  # the derivatives wrt to the other parameters.
  # It returns a objective function to be passed to numDeriv::jacobian
  derWrapCreator <- function(ii, kp, pr, derFun){
    outFun <- function(x){
      if( length(kp) < length(pr) ) { x <- c(pr[-kp], x) }
      derFun(x)[[ii]]
    }
    return(outFun)
  }
  
  d <- length( param )
  out <- list()
  
  for(od in ord){
    # In the following:
    # keep is a list, where each entry is a vector of indexes of parameters wrt which we want to differentiate
    # tot is the total number of derivatives (d first gradient, (d^2 + 1)/2 for Hessian and so on...)
    
    keep <- list()
    
    if( od == 1 ){
      fd <- drop( jacobian(func = obj$d0, x = param) )
    } else {
      
      if(od == 2){
        derFun <- obj$d1
        tot <- d
        for(ii in 1:d){ keep[[ii]] <- ii:d }
      }
      
      if(od == 3){
        derFun <- obj$d2
        tot <- (d^2 + d)/2
        kk <- 1
        keep <- list()
        for(ii in 1:d){
          for(jj in ii:d){
            keep[[kk]] <- jj:d
            kk <- kk + 1
          }
        }
      }
      
      fd <- list()
      for(ii in 1:tot){ 
        wr <- derWrapCreator(ii = ii, kp = keep[[ii]], pr = param, derFun = derFun)
        fd[[ii]] <- jacobian(func = wr, x = param[keep[[ii]]])
      }
      fd <- unlist(fd)
      
    }
    
    exd <- drop( obj[[paste0("d", od)]](param) )
    derPair <- cbind(exd, fd)
    colnames(derPair) <- c("EX", "FD")
    
    out[[paste0("fd", od)]] <- matrix(unlist(derPair), ncol = ncol(derPair))
    
  }
  
  return(out)
  
}