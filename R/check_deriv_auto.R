#'
#' Checking derivatives by finited differences 
#' 
#' @description Function for comparing exact and numerical derivatives.
#' @param np number of parameter vectors to be simulated.
#' @param simfun function to simulate parameter vectors.
#' @param obj for instance the output of createSH
#' @param ord the order to derivatives to be checked.
#' @param trans transformation to be applied to the derivatives
#' @param ... arguments passed to \code{obj$initialize}.
#' @param fun the name of the slot of \code{obj} containing the derivative function.
#' @name check_deriv_auto
#' @rdname check_deriv_auto
#' @export check_deriv_auto
#'
check_deriv_auto <- function(np, simfun, obj, ord = 1, trans = identity, ..., fun = "derobj"){
  
  pars <- simfun(np)
  
  # Simulate observations, recreate object and evaluate derivatives
  der <- list()
  for(ii in 1:np){
    obj <- obj$initialize(param = pars[ii, , drop = TRUE], ...)
    der[[ii]] <- check_deriv(obj = wrap_der_fun(obj[[fun]], dropH = TRUE), param = pars[ii, ], ord = ord)
  }
  
  # Gather all simulations with same derivative order
  out <- list()
  for(jj in ord){
    out[[paste0("d", jj)]] <- do.call("rbind", lapply(der, "[[", paste0("fd", jj)))
  }
  
  for(kk in 1:length(out)){
    plot(trans(out[[kk]][ , 1]), trans(out[[kk]][ , 2]), main = names(out)[[kk]], 
         xlab = colnames(out[[1]])[1], ylab = colnames(out[[1]])[2])
    abline(0, 1, col = 2)
  }
  
  out$param <- pars
  
  return(out)
}
