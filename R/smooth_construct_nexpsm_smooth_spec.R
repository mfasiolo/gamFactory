#'
#' Nested adaptive exponential smoothing effects for mgcv
#' 
#' @name smooth.construct.nexpsm.smooth.spec
#' @rdname smooth.construct.nexpsm.smooth.spec
#' @importFrom Matrix rankMatrix
#' @export
#'
smooth.construct.nexpsm.smooth.spec <- function(object, data, knots)
{ 
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # Extract variables from data frame:
  # - "y" is variable to be exponentially smoothed (stored in vector or matrix)
  # - "x" is model matrix used to model the exp smoothing rate
  # - "times" is the vector of times that at which the response variable of the 
  #           GAM is observed (hence not the same at the GAM above). 
  Xi <- data[[object$term]]
  n <- nrow( Xi )
  nms <- colnames(Xi)
  x <- as.vector( Xi[ , which(nms == "y")] )
  times <- NULL
  tmp <- which(nms == "times")
  if( length(tmp) ){
   times <- Xi[ , tmp]
  }
  Xi <- Xi[ , which(nms == "x"), drop = FALSE]
  nrep <- ceiling( length(x)/n )
  Xi <- lapply(0:(nrep-1), function(ii){
    dXi <- ncol(Xi) / nrep
    as.matrix(Xi[ , (ii*dXi + 1):(dXi*(ii+1))])
  })
  Xi <- do.call("rbind", Xi)
  Si <- si$S
  di <- ncol(Xi) + 1 # + 1 because we have extra scaling parameter multiplying output of the exp smooth.
                     # Note that this extra parameter should NOT be penalised by Si
  
  if( !is.null(Si) ){ # Reparametrise Xi so that the penalty on alpha is diagonal
    si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = rankMatrix(Si)))
  } else { # No inner penalty 
    si$X <- Xi
    si$B <- diag(nrow = ncol(Xi))
    si$rank <- 0 
  }
  
  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = 1 (target variance)
    g <- expsmooth(y = x, Xi = si$X, beta = rep(0, di-1), times = times)$d0
    alpha <- si$alpha <- c(log(1/sd(g)), rep(0, di-1))
  } else {
    alpha <- solve(si$B) %*% alpha
    g <- expsmooth(y = x, Xi = si$X, beta = alpha, times = times)$d0
    alpha <- si$alpha <- c(log(1/sd(g)), alpha)
  }
  
  # Center and scale the initialized inner linear predictor
  data[[object$term]] <- exp(alpha[1]) * (g - mean(g))
  
  si$x <- x
  si$times <- times
  
  out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # Add inner penalty matrix: rbind(0, cbind(0, si$S) makes so first element of alpha (inner parameters)
  # is unpenalised (this is the scaling parameter)
  if( !is.null(Si) ){
    dsmo <- out$bs.dim - di
    si <- out$xt$si
    out$S[[2]] <- rbind(cbind(rbind(0, cbind(0, si$S)), matrix(0, di, dsmo)),
                        cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo)))
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank)
    out$rank <- c(out$rank, si$rank)
  }
  
  class(out) <- c("nexpsm", "nested")
  return( out )
} 
