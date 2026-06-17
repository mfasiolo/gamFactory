#'
#' Compute jacobian of nested effects
#' 
#' @noRd
get_jacobian.nested <- function(object,data, param){
  
  if(class(object)[1] == "si"){
    return( .get.jacobian.si(object, data, param) ) 
  }
  if(class(object)[1] == "nexpsm" || class(object)[1] == "mgks"){
    return( .get.jacobian.nexpsm.mgks(object, data, param) )
  }
  if(class(object)[1] == "si_nexpsm"){ 
    return(.get.jacobian.si_nexpsm(object, data, param))
  }
  
  stop("I do not know this effect type")
  
}


# ------------------------------subfunction------------------------------
.get.jacobian.si <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  
  JJ <- cbind(drop(store$X1 %*% beta) * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma
  return(list("JJ" = JJ, "xa" = NULL) )
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
  
  return(list("JJ" = JJ, "xa" = NULL) )
}

.get.jacobian.si_nexpsm <- function(object, data, param){
  na <- length(object$xt$si$alpha)
  beta <- param[ -(1:na) ]
  
  x_nest <- gamFactory:::Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  X1beta <- drop(store$X1 %*% beta)  # derivative to beta
  
  # ∂eta/∂alpha_si = (dM/dz %*% beta) * ∂z/∂alpha_si
  J_alpha_si <- X1beta * x_nest$xa_dalpha_si  # n × na_si
  # ∂eta/∂alpha_nexp = (dM/dz %*% beta) * ∂z/∂alpha_nexp
  J_alpha_nexp <- X1beta * x_nest$xa_dalpha_nexp  # n × na_nexp
  # ∂eta/∂beta = M
  J_beta <- store$X0
  
  JJ <- cbind(J_alpha_nexp, J_alpha_si, J_beta)
  
  return( list("JJ" = JJ, 
               "z" = x_nest$z,
               "xa" = x_nest$xa, 
               "xsm_raw" = x_nest$xsm_raw) )
}









