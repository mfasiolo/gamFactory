#'
#' Predict using nested effects
#' 
#' @name Predict.matrix.nested
#' @rdname Predict.matrix.nested
#' @export
#'
#'
Predict.matrix.nested <- function(object, data, ...){
  
  if(class(object)[1] == "si"){
    return( .predict.matrix.si(object, data, ...) ) 
  }
  if(class(object)[1] == "nexpsm"){
    return( .predict.matrix.nexpsm(object, data, ...) ) 
  }
  if(class(object)[1] == "mgks"){
    return( .predict.matrix.mgks(object, data, ...) ) 
  }
  if(class(object)[1] == "si_nexpsm"){
    return( .predict.matrix.si_nexpsm(object, data, ...) )  # <-new type
  }
  
  stop("Predict.matrix.nested --- I do not know this effect type")
  
}