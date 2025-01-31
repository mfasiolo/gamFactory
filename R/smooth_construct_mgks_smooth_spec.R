#'
#' Nested MGKS effects for mgcv
#' 
#' @name smooth.construct.mgks.smooth.spec
#' @rdname smooth.construct.mgks.smooth.spec
#' @importFrom mgcv smooth.construct
#' @examples 
#' n <- 100
#' n0 <- 50
#' p <- 5
#' dat <- data.frame(y = rnorm(n))
#' dat$Xi <- matrix(rnorm(n*p), n, p)
#' X0 <- matrix(rnorm(n0*p), n0, p)
#' x <- rnorm(n0)
#' aaa <- smoothCon(s(Xi, bs = "mgks", m = c(3, 2),
#'             xt = list(si = list(X0 = X0, x = x), sumConv = FALSE)),
#'           data = dat)
#' @export
smooth.construct.mgks.smooth.spec <- function(object, data, knots)
{ 
  si <- object$xt$si
  
  # We have n0 original locations (one on each row of X0) and corresponding observations in y0
  Xi <- data[[object$term]]
  n <- nrow( Xi )
  nms <- colnames(Xi)
  y0 <- Xi[ , which(nms == "y"), drop = FALSE]
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
  d <- kk
  
  si$x <- y0
  si$dist <- Dist
  
  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = 1 (target variance)
    # Other elements of alpha set to the negative marginal standard deviations / 10.
    # Dividing by 10 seems a good compromise between under- and over-smoothing.
    g <- mgks(y = si$x, dist = Dist, beta = -log(sapply(si$dist, sd)/10))$d0
    alpha <- si$alpha <- c(log(1/sd(g)), -log(sapply(si$dist, sd)/10)) 
  } else {
    g <- mgks(y = si$x, dist = Dist, beta = alpha)$d0
    alpha <- si$alpha <- c(log(1/sd(g)), alpha)
  }
  
  # Center and scale the initialized inner linear preditor
  data[[object$term]] <- exp(si$alpha[1]) * (g - mean(g))
  
  out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  class(out) <- c("mgks", "nested")
  return( out )
} 

