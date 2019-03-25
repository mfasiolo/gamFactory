#'
#' Checking derivatives by finited differences 
#' 
#' @description Function for comparing exact and numerical derivatives.
#' @param np number of parameter vectors to be simulated.
#' @param parSim function to simulate parameter vectors.
#' @param obj for instance the output of createSH
#' @param ord the order to derivatives to be checked.
#' @param trans trasformation to be applied to the derivatives
#' @param ... arguments passed to \code{obj$initialize}.
#' @param derFunNam the name of the slot of \code{obj} containing the derivative function.
#' @name derivCheck
#' @rdname derivCheck
#' @export derivCheck
#'
derivCheck <- function(np, parSim, obj, ord = 1, trans = identity, ..., derFunNam = "derObj"){
  
  pars <- parSim(np)
  
  # Simulate observations, recreate object and evaluate derivatives
  der <- list()
  for(ii in 1:np){
    obj <- obj$initialize(param = pars[ii, , drop = TRUE], ... )
    der[[ii]] <- fdDeriv(obj = derFunWrapper(obj[[derFunNam]], dropHessian = TRUE), param = pars[ii, ], ord = ord)
  }
  
  # Gather all simulations with same derivative order
  out <- list()
  for(jj in ord){
    out[[paste0("d", jj)]] <- do.call("rbind", lapply(der, "[[", paste0("fd", jj)))
  }
  
  par(mfrow = c(2, 2))
  for(kk in 1:length(out)){
    plot(trans(out[[kk]][ , 1]), trans(out[[kk]][ , 2]), main = names(out)[[kk]], 
         xlab = colnames(out[[1]])[1], ylab = colnames(out[[1]])[2])
    abline(0, 1, col = 2)
  }
  
  return(out)
}
