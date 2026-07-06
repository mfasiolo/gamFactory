#'
#' Comparing analytic with numerical finite differences
#'
#' @description Function for getting higher derivatives by finite differences.
#'              Any order is supported (1st order derivatives are checked against
#'              \code{obj$d0}, 2nd order against \code{obj$d1}, and so on up to
#'              \code{obj$d\{max(ord)\}} against \code{obj$d\{max(ord) - 1\}}), as long as
#'              \code{obj} provides the corresponding \code{d0, ..., d\{max(ord)\}} entries,
#'              each following the \code{[gamFactory::llk_gaussian]}-style convention of
#'              returning a list whose elements are ordered by non-decreasing parameter
#'              index tuples (e.g. \code{d2} ordered as (1,1), (1,2), ..., (1,d), (2,2), ...).
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

  # Generate all non-decreasing tuples of length `len` with entries in 1:d, in
  # the same nested-loop order (i_1 in 1:d, i_2 in i_1:d, ...) used throughout
  # gamFactory to lay out the elements of d1, d2, d3, d4, ... (e.g. mgcv::trind.generator).
  nonDecreasingTuples <- function(d, len){
    if( len == 0 ){ return( list(integer(0)) ) }
    tuples <- list()
    rec <- function(prefix, start){
      if( length(prefix) == len ){
        tuples[[length(tuples) + 1]] <<- prefix
      } else {
        for(v in start:d){ rec(c(prefix, v), v) }
      }
    }
    rec(integer(0), 1)
    return(tuples)
  }

  d <- length( param )
  out <- list()

  for(od in ord){
    # In the following:
    # keep is a list, where each entry is a vector of indexes of parameters wrt which we want to differentiate

    if( od == 1 ){
      fd <- drop( jacobian(func = obj$d0, x = param) )
    } else {

      derFun <- obj[[paste0("d", od - 1)]]

      # Tuples of length (od - 1) indexing the elements of derFun's output, e.g.
      # single indices for od = 2 (derFun = d1), pairs i<=j for od = 3 (derFun = d2),
      # triples i<=j<=k for od = 4 (derFun = d3), and so on.
      tuples <- nonDecreasingTuples(d, od - 1)

      keep <- lapply(tuples, function(.tup) .tup[length(.tup)]:d)

      fd <- list()
      for(ii in seq_along(tuples)){
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