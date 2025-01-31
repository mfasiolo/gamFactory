#'
#' Compute jacobian of nested effects
#' 
#' @noRd
.get_jacobian.nested <- function(object,data, param){
  
  if(class(object)[1] == "si"){
    return( .get.jacobian.si(object, data, param) ) 
  }
  if(class(object)[1] == "nexpsm" || class(object)[1] == "mgks"){
    return( .get.jacobian.nexpsm.mgks(object, data, param) )
  }
  
  stop("I do not know this effect type")
  
}

.get.jacobian.si <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  
  JJ <- cbind(drop(store$X1 %*% beta) * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma
  return(JJ)
}

.get.jacobian.nexpsm.mgks <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  X1beta <- drop(store$X1 %*% beta)
  
  JJ <- cbind(X1beta * x_nest$xa, # df/da = M1%*%b * ds/da0 (where ds/da0 = s because s = exp(a0) xa)
              X1beta * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma 
  
  return(JJ)
}











