#' Predict fitted values at all time points from a nested GAM-NL model
#'
#' smoothing components.
#'
#' @param object A fitted \code{gamnl} model object obtained from
#'        \code{\link[gamFactory]{gamnl}} or similar fitting functions.
#' @param X_ori A matrix of covariates corresponding to the inner
#'        linear transformation component.
#' @param w_ori A matrix of covariates corresponding to the outer
#'        exponential smoothing component.
#'
#' @return
#' A numeric vector containing the predicted values for all time points
#' in the same order as the input rows of \code{X_ori} and \code{w_ori}.
#'
#'
#' @export

pre_all_times <- function(object,
                          X_ori,
                          w_ori){
  newdata <- data.frame(
    xw = I(cbind(X_ori,w_ori))
  )
  
  pre_all <- gamFactory:::predict.gamnl(
    object      = object,
    newdata    = newdata
  )
  
  pre_all <- pre_all[,1]
  
  return(pre_all)
  
}