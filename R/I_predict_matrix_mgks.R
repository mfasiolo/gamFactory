#'
#' Predict using MGKS smooth effects
#' 
#' @noRd
.predict.matrix.mgks <- function(object, data, get.xa = FALSE){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  # Do NOT remove this! We need it in postproc_gam_nl.
  if( is.null(si$xm) ){
    si$xm <- 0
  }
  
  # We have n0 original locations and corresponding observations in y0
  Xi <- data[[object$term]]
  n <- nrow( Xi )
  nms <- colnames(Xi)
  y0 <- Xi[ , which(nms == "y")]
  if( !ncol(y0) ){
    y0 <- si$y0
  }
  Dist <- list()
  kk <- 1
  while( TRUE ){
    idx <- which(startsWith(nms, "d") & endsWith(nms, as.character(kk)) & sapply(nms, function(.x) nchar(.x) == 2))
    if( !length(idx) ){
      break
    }
    Dist[[kk]] <- Xi[ , idx, drop = FALSE]
    kk <- kk + 1
  }
  
  alpha <- si$alpha
  a0 <- alpha[1]
  a1 <- alpha[-1]
  
  xsm_list <- mgks(y = y0, dist = Dist, beta = a1, deriv = get.xa)
  xsm_unscaled <- xsm_list$d0 - si$xm
  xsm <- exp(a0) * xsm_unscaled 
  
  if(get.xa){
    return(list(xa = xsm, xa_da = exp(a0)*xsm_list$d1))
  }
  
  # Compute outer model matrix
  X1 <- object$xt$basis$evalX(x = xsm, deriv = 0)$X0
  
  # Total model matrix is X1 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X1), length(alpha)), X1) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xsm_unscaled
  
  return(Xtot)
  
}