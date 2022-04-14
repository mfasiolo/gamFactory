#'
#' Nested MGKS effects for mgcv
#' 
#' @name smooth.construct.mgks.smooth.spec
#' @rdname smooth.construct.mgks.smooth.spec
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
  
  Xi <- data[[object$term]]
  di <- (ncol(Xi)-1)/2 + 1
  n <- nrow( Xi )
  n0 <- si$n0
  si$x <- x <- Xi[1:n0, 1]
  si$X0 <- X0 <- Xi[1:n0, 2:di, drop = FALSE]
  si$X <- Xi <- Xi[ , -(1:di), drop = FALSE]
  
  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = 1 (target variance)
    # Other elements of alpha set to the negative marginal standard deviations / 10.
    # Dividing by 10 seems a good compromise between under- and over-smoothing.
    g <- mgks(y = x, X = Xi, X0 = X0, beta = -log(colSds(X0)/10))$d0
    alpha <- si$alpha <- c(log(1/sd(g)), -log(colSds(X0)/10)) 
  } else {
    g <- mgks(y = x, X = Xi, X0 = X0, beta = alpha[-1])$d0
  }
  
  # Center and scale the initialiized inner linear preditor
  data[[object$term]] <- exp(alpha[1]) * (g - mean(g))
  
  out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)

  class(out) <- "mgks.smooth"
  return( out )
} 
