#'
#' Defining nested smooths in GAM formulae
#' 
#' @name s_nest
#' 
#' @description Function used to defined nested smooth effects. It works similarly to 
#'              \link{mgcv::s}.
#'              
#' @param ... a list of variables that are used to defined the nested smooth effect.
#' @param trans the type of covariate transformation to use (see, e.g., \link{trans_mgks}).
#' @param k the number of basis functions used to construct the outer smooth effect. 
#' @param m a vector of two positive integers. The first indicates the order of the 
#'          outer B-spline basis (see \link{basis_bspline}), the second is the order of the corresponding P-spline
#'          penalty. 
#' @param sarg a list of further arguments to be passed to the effect constructor. 
#' @rdname s_nest
#' @export
#'
s_nest <- function(..., trans, k, m, sarg){
  
  vars <- as.list(substitute(list(...)))[-1]
  
  if( missing(sarg) ){ sarg <- list() }

  if( is.null(sarg$xt) ){ sarg$xt <- list() }
  sarg$xt$sumConv <- FALSE
  sarg$bs <- trans$type
  trans$type <- NULL
  if( length(trans) ){ sarg$xt$si <- trans }
  
  if( missing(k) ){
    k <- 10 + 7 # Because we impose 7 constraints
  }
  sarg$k <- k 
  
  if( missing(m) ){ m <- c(4, 2) }
  if( length(m) == 1 ){ m <- c(m, 2) }
  if( m[1] < 4 ){ message("m[1] should be >= 4, lower values might lead to convergence issues.") }
  sarg$m <- m
  
  # Call mgcv::s() to create smooth
  sm <- do.call("s", c(vars, sarg))
  
  # Create formula for term "s(...)" that can be used to build GAM model formula
  vnam <- as.character(vars)
  form <- do.call("call", c("s", c(vnam, sarg)))
  form <- paste(deparse(form), collapse="")
  
  # Transform "s(\"var1\", \"var2\", ...)" into "s(var1, var2, ...)" 
  for(ii in 1:length(vnam)){
    form <- gsub("\\s+", " ", gsub(paste0("\"",vnam[ii],"\""), vnam[ii], form, fixed=TRUE))
  }
  
  sm$form_term <- form
  
  return(sm)
  
}