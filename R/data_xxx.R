#' Prepare data matrix for a nested smooth effect with \code{trans = trans_linear}
#'
#' Builds the matrix to be passed as a column of the data frame used to fit a \code{gam_nl}
#' model with a single index smooth term, i.e. a smooth term 
#' \eqn{s(\tilde{s}(\mathbf{x}_i)) = s(\mathbf{x}_i^\top \mathbf{a})}.
#'
#' @usage data_linear(data, cnames = NULL)
#'
#' @param data List or matrix containing the variables to be linearly combined.
#'   If a list, each element must be a numeric vector of the same length; these
#'   are column-bound to form the output matrix.  If already a matrix, it is
#'   returned as-is (with column names optionally replaced by \code{cnames}).
#' @param cnames Optional character vector of length \code{p} giving the column
#'   names of the output matrix, where \code{p} is the number of variables
#'   (list elements or matrix columns) in \code{data}.  If \code{NULL}
#'   (default), existing column names are preserved, or columns are left
#'   unnamed if none are present.
#'
#' @return A numeric matrix with \code{n} rows and \code{p} columns, where
#'   \code{n} is the length of each variable and \code{p} is the number of
#'   variables in \code{data}.  Column names are set to \code{cnames} when
#'   supplied.
#'
#' @examples
#' set.seed(6436)
#' n <- 1000; p <- 3; b <- 1:p
#' X <- lapply(1:p, \(ii) rnorm(n)) # equivalently matrix(runif(p * n), ncol = p)
#' cn <- paste0("x.", 1:p)
#'
#' Xb <- do.call("cbind", X) %*% b
#' y <- 2 * sin(Xb) + rnorm(n)
#' dat <- data.frame(y = y)
#'
#' dat$X <- data_linear(X, cnames = cn)
#' str(dat)
#' colnames(dat$X)
#'
#' @export
data_linear <- function(data, cnames = NULL){
  res <- if(!is.matrix(data)) do.call("cbind", data) else data
  if(!is.null(cnames)) colnames(res) <- cnames
  res
}

#' Prepare data matrix for a nested smooth effect with \code{trans = trans_nexpsm}
#'
#' Builds the matrix to be passed as a column of the data frame used to fit a \code{gam_nl}
#' model with a nested exponential smoothing term, i.e. a smooth term 
#' \eqn{\tilde{s}(x_i) = \tilde{s}_i = \omega_{i}\tilde{s}_{i-1}+(1-\omega_{i})x_{i}}, 
#' for \eqn{i \geq 1,\ \tilde{s}_0 = x_0,\ \omega_i = \phi(\tilde{\mathbf{x}}_i^\top {\bm a})}, 
#' where \eqn{\phi} is the logistic function, \eqn{\tilde{\mathbf{x}}_i} is a given vector.
#' 
#' @usage data_nexpsm(y, X=NULL, times = NULL, n_obs=NULL)
#'
#' @param y Numeric vector of length \code{m} or matrix with \code{m} rows.
#'   The series to be exponentially smoothed.  it can be of the same length as 
#'   the number of observations in the response variable, or longer than it. 
#'   If \code{y} is a matrix, it must have a single column.  
#' @param X Numeric vector or matrix with \code{m} rows (same \code{m} as
#'   \code{y}).  Model matrix for the exponential smoothing rate \eqn{\alpha}.
#'   Include an intercept column of ones explicitly if a constant rate is desired.  
#'   Defaults to \code{NULL}, in which case a single intercept column is used.
#' @param times Optional numeric vector of length \code{n_obs} giving the
#'   observation times of the response variable.  When supplied, it is
#'   appended as a column named \code{"times"}. 
#'   If \code{NULL} (default) the smoother assumes regular unit spacing.
#' @param n_obs Optional integer giving the number of response
#'   observations.  Used to determine \code{nrep = ceiling(nrow(y) / n_obs)}.
#'   If both \code{n_obs} and \code{times} are \code{NULL}, defaults to
#'   \code{nrow(y)} (i.e. \code{nrep = 1}).  If \code{times} is supplied,
#'   \code{n_obs} is set to \code{length(times)} and this argument is ignored.
#' @return A numeric matrix with \code{n_obs} rows and
#'   \code{nrep + nrep * ncol(X) + as.integer(!is.null(times))} columns,
#'   with column names set to the \code{"y"} / \code{"x"} / \code{"times"}
#'   convention required by \code{smooth.construct.nexpsm.smooth.spec}.
#'   Padding rows (if \code{nrow(y)} is not an exact multiple of \code{n_obs})
#'   are filled with \code{Inf}.
#'
#' @examples
#' ### Intercept only, regular spacing
#' library(gamFactory)
#' n   <- 200
#' Xi <- matrix(1, ncol = 1, nrow=n)
#' x <- seq(0, 2 * pi, length.out = n)
#' wi <- sin(x) + rnorm(length(x), 0, 0.5)
#' z <- expsmooth(y = wi, Xi = Xi, beta = c(-0.8))$d0
#' dat <- data.frame(y = z + rnorm(n))
#' dat$mat <- data_nexpsm(y=wi)
#' str(dat)
#' head(dat)
#'
#' ### With observation times
#' pp <- 1.5
#' times <- cumsum(rpois(n, lambda = 2))
#' Xi <- matrix(1, ncol = 1, nrow=pp*n)
#' x <- seq(0, 2*pi, length.out = pp*n)
#' wi <- sin(x) + rnorm(length(x), 0, 0.5)
#' z <- expsmooth(y = wi, Xi = Xi, beta = c(-0.8))$d0[times]
#' dat <- data.frame(y = z + rnorm(n))
#' dat$mat <-  data_nexpsm(y=wi, times=times)
#' str(dat)
#' head(dat)
#' 
#' ### With observation times and linear predictors for the smoothing rate
#' X <- cbind(1, matrix(rnorm(2*length(x)), ncol = 2))
#' z <- expsmooth(y = wi, Xi = X, beta = c(-0.8, 1, 0.4))$d0[times]
#' dat <- data.frame(y = z + rnorm(n))
#' dat$mat <- data_nexpsm(y=wi,X=X, times=times)
#' str(dat)
#' head(dat)
#' @export
data_nexpsm <- function(y, X=NULL, times = NULL, n_obs=NULL) {
    
  # variable to be smoothed
  if (is.vector(y)) {y <- matrix(y, ncol = 1)} else if (!is.matrix(y)) {stop("`y` must be a numeric vector or matrix.")
    } else if (ncol(y) > 1) {stop("`y` must be a vector or a matrix with a single column.")}
  ny    <- nrow(y)
    
  # variables for the inner linear predictor of the smoothing rate
  # if X is null, then exponential smoothing with a constant rate is fitted.  
  # If X is a vector, it is coerced to a single-column matrix.  
  # If X is a matrix, it is used as-is.  In all cases X must have the same number of rows as y.
  if (is.null(X)) {
    X <- matrix(1, nrow = ny, ncol = 1)
  } else if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  } else if (!is.matrix(X)) {
    stop("`X` must be a numeric vector or matrix.")
  }
  if (nrow(X) != ny) {stop("`X` must have the same number of rows as `y`.","\n  nrow(y) = ",ny,", nrow(X) = ",nrow(X))}

  if (!is.null(times)) {
    n_obs <- length(times)
  } else if (is.null(n_obs)) {
    n_obs <- ny
  }
  nrep <- ceiling(ny/n_obs)
  
  # Pad y and X to nrep * n_obs rows with Inf
  pad <- nrep * n_obs - ny
  y <- matrix(c(y, rep(Inf, pad * ncol(y))), ncol = ncol(y))
  X <- rbind(X, matrix(Inf, nrow = pad, ncol = ncol(X)))
  
  # Reshape into (n_obs x nrep) and (n_obs x nrep*ncol(X))
  # Both y and X are filled row-by-row (byrow = TRUE).  The constructor later
  # unpacks with as.vector(t(Xi[, cols])), which reads column-by-column from
  # the transposed matrix.
  mat <- list(matrix(as.vector(y),n_obs, nrep,byrow = TRUE),
              matrix(as.vector(t(X)),n_obs, nrep*ncol(X),byrow = TRUE)
  )
  cnames <- c(rep("y", nrep), rep("x", nrep * ncol(X)))
  
  # Optionally append times
  if (!is.null(times)) {
    mat[[3]] <- times
    cnames   <- c(cnames, "times")
  }
  
  res <- do.call("cbind", mat)
  colnames(res) <- cnames
  res
}


#' Prepare data matrix for a double nested smooth effect with
#' \code{trans = trans_linear_nexpsm}
#'
#' Builds the matrix to be passed as a column of the data frame used to fit a
#' \code{gam_nl} model with a double nested effect. The inner layer first
#' constructs a single-index variable
#' \eqn{z_i = \mathbf{x}*i^\top {\bm a}*{si}}, and the outer layer then applies
#' exponential smoothing with a possibly covariate-dependent smoothing rate,
#' \eqn{\tilde{s}*i = \omega_i \tilde{s}*{i-1} + (1-\omega_i) z_i},
#' for \eqn{i \geq 1}, with initial value \eqn{\tilde{s}_0 = z_0}. The smoothing
#' weight is given by
#' \eqn{\omega_i = \phi(\mathbf{w}*i^\top {\bm a}*{nexp})}, where \eqn{\phi}
#' is the logistic function. Here \eqn{\mathbf{x}_i} is the covariate vector for
#' the single-index layer and \eqn{\mathbf{w}_i} is the covariate vector for the
#' exponential-smoothing-rate layer.
#'
#' @usage data_linear_nexpsm(X_si, X_nexp = NULL, times = NULL,
#'   n_obs = NULL)
#'
#' @param X_si Numeric matrix with \code{n} rows and more than one column.
#'   The covariates used in the inner single-index transformation
#'   \eqn{z_i = \mathbf{x}*i^\top {\bm a}*{si}}. Each row corresponds to one
#'   time point in the underlying full sequence.
#' @param X_nexp Optional numeric vector or matrix. The covariates used in the
#'   linear predictor of the exponential smoothing rate
#'   \eqn{\omega_i = \phi(\mathbf{w}*i^\top {\bm a}*{nexp})}. If a vector is
#'   supplied, it is coerced to a single-column matrix. If \code{NULL}, a single
#'   column of ones is used, corresponding to a constant smoothing rate.
#' @param times Optional numeric vector giving the observation times of the
#'   response variable. When supplied, \code{n_obs} is set to
#'   \code{length(times)} and the supplied \code{n_obs} is ignored. Unlike
#'   \code{data_nexpsm()}, the \code{times} vector is not appended to the
#'   returned matrix. It is used only to determine the number of response
#'   observations, and, when \code{time_gap = TRUE}, to construct a time-gap
#'   covariate for the exponential-smoothing-rate layer.
#' @param n_obs Optional integer giving the number of response observations.
#'   Used to determine \code{nrep = ceiling(nrow(X_si) / n_obs)}. If both
#'   \code{n_obs} and \code{times} are \code{NULL}, defaults to
#'   \code{nrow(X_si)}, i.e. \code{nrep = 1}. If \code{times} is supplied,
#'   this argument is ignored.
#'
#' @return A numeric matrix with \code{n_obs} rows. The first block of columns,
#'   named \code{"x"}, contains the reshaped single-index covariates
#'   \code{X_si}. The second block of columns, named \code{"w"}, contains the
#'   reshaped exponential-smoothing-rate covariates \code{X_nexp}. Padding rows, 
#'   used when \code{nrow(X_si)} is
#'   not an exact multiple of \code{n_obs}, are filled with \code{Inf}.
#'
#' @examples
#' ### Covariate-dependent smoothing rate, regular spacing
#' library(gamFactory)
#' set.seed(6436)
#'
#' n <- 1000
#' n_si <- 3
#' n_nexp <- 2
#'
#' X_si   <- matrix(rnorm(n * n_si), n, n_si)
#' X_nexp <- matrix(rnorm(n * n_nexp), n, n_nexp)
#'
#' ## Parameter order: [alpha_nexp, alpha_si]
#' param <- rnorm(n_nexp + n_si)
#'
#' y <- gamFactory:::deriv_si_nexp(
#'   X_si = X_si,
#'   X_nexp = X_nexp,
#'   param = param,
#'   deriv = 0
#' )$d0
#'
#' dat <- data.frame(y = y + rnorm(n))
#' dat$X <- data_linear_nexpsm(
#'   X_si = X_si,
#'   X_nexp = X_nexp
#' )
#' str(dat)
#' head(dat)
#'
#' ### With observation times
#' time_seq <- sort(sample(seq_len(n), 500))
#'
#' dat <- data.frame(y = y[time_seq] + rnorm(length(time_seq)))
#' dat$X <- data_linear_nexpsm(
#'   X_si = X_si,
#'   X_nexp = X_nexp,
#'   times = time_seq
#' )
#' str(dat)
#' head(dat)
#'
#' ### Constant smoothing rate
#' dat <- data.frame(y = y + rnorm(n))
#' dat$X <- data_linear_nexpsm(
#'   X_si = X_si,
#'   X_nexp = NULL
#' )
#' str(dat)
#' head(dat)
#'
#'
#' @export
data_linear_nexpsm <- function(X_si, X_nexp = NULL, 
                               times = NULL, n_obs = NULL){
  if (!is.matrix(X_si) || ncol(X_si) == 1) stop("`X_si` must be a matrix with more than 1 column.")
  n <- nrow(X_si)
  
  if (!is.null(times)) {
    if(!is.null(n_obs) && n_obs!= length(times)) message("`n_obs` is ignored when `times` is supplied.  Using `n_obs = length(times)`.")
    n_obs <- length(times)
  } else if (is.null(n_obs)) {
    n_obs <- n
  }
  
  if(!is.null(X_nexp)){
    if (is.vector(X_nexp)) { X_nexp <- matrix(X_nexp, ncol = 1)
    } else if (!is.matrix(X_nexp)) stop("`X_nexp` must be a numeric vector or matrix.")
  } else {X_nexp <- matrix(1, nrow = n, ncol = 1) }
  
  nrep <- ceiling(n / n_obs)
  # Pad to nrep * n_obs rows with Inf 
  pad    <- nrep * n_obs - n
  X_si   <- rbind(X_si, matrix(Inf, pad, ncol(X_si)))
  X_nexp <- rbind(X_nexp, matrix(Inf, pad, ncol(X_nexp)))
  
  n_si <- ncol(X_si)
  n_nexp <- ncol(X_nexp)

  # Reshape into (n_obs x nrep*n_si) and (n_obs x nrep*n_nexp)
  # Filled row-by-row (byrow = TRUE); constructor unpacks with
  # as.vector(t(Xi[, cols])) which reads column-by-column from the transpose.
  si_block   <- matrix(as.vector(t(X_si)),n_obs, nrep * n_si, byrow = TRUE)
  nexp_block <- matrix(as.vector(t(X_nexp)), n_obs, nrep * n_nexp, byrow = TRUE)
  
  # Assemble final matrix
  cnames <- c(rep("x", ncol(si_block)), rep("w", ncol(nexp_block)))
  mat    <- list(si_block, nexp_block)
  
  res <- do.call("cbind", mat)
  colnames(res) <- cnames
  res
}


#' Prepare data matrix for a nested smooth effect with \code{trans = trans_mgks}
#'
#' Builds the matrix to be passed as a column of the data frame used to fit a \code{gam_nl}
#' model with a nested Multivariate Gaussian Kernel Smoothing (MGKS) term, i.e. a smooth term
#' \deqn{\tilde{s}({\bf x}_i) = \frac{\sum_{j\in \mathcal{N}_i}K_{{\bm a}}({\bf x}_{i}, 
#' {{\bf x}}_{j})z_{j}}{\sum_{q\in \mathcal{N}_i}K_{{\bm a}}({\bf x}_{i},{{\bf x}}_{q})},}
#' where \eqn{K_{\mathbf{a}}} is a kernel function parameterized by \eqn{\mathbf{a}},
#' \eqn{z_{ij}} are scalar observations at locations \eqn{\mathbf{x}_{ij}}, and
#' \eqn{\mathbf{x}_i} is the target location at which the response is observed.
#'
#' @usage data_mgks(y, dist)
#'
#' @param y Numeric matrix with \code{n} rows and \code{p} columns, where \code{n} is
#'   the number of response observations and \code{p} is the number of neighbouring
#'   locations. Each column contains the scalar covariate \eqn{z_{ij}} observed at
#'   the \eqn{j}-th neighbouring location.
#' @param dist List of numeric matrices, each with \code{n} rows and \code{p} columns.
#'   Each matrix contains the pairwise distances (or squared distances) along one
#'   spatial dimension between the \code{n} target locations and the \code{p}
#'   neighbouring locations. The length of the list determines the number of spatial
#'   dimensions \eqn{d}. All matrices must have \code{p} columns, matching \code{ncol(y)}.
#'
#' @return A numeric matrix with \code{n} rows and \code{p * (1 + d)} columns, where
#'   \code{d = length(dist)}, with column names set to the \code{"y"} / \code{"d1"} /
#'   \code{"d2"} / \ldots convention required by
#'   \code{smooth.construct.mgks.smooth.spec}: \code{p} columns named \code{"y"}
#'   followed by \code{p} columns named \code{"d1"}, \code{p} columns named
#'   \code{"d2"}, and so on.
#'
#' @examples
#' ### Two-dimensional kernel smoothing over n0 neighbouring locations
#' set.seed(123)
#' n0 <- 50; n <- 1000
#' # locations where z is observed
#' X0 <- cbind(runif(n0, -1, 1), runif(n0, -4, 4))
#' trueF <- function(x) 3 * x[ , 1] + x[ , 2]^2
#' z <- t(sapply(1:n, function(ii) trueF(X0)))
#'
#' # locations where response variable and the other covariates are observed
#' X <- cbind(runif(n, -1, 1), runif(n, -4, 4))
#' 
#' # compute distance between each target location and each neighbouring location along each spatial dimension
#' dist <- lapply(1:ncol(X), function(dd){
#'   t(sapply(1:nrow(X), function(ii) (X[ii, dd] - X0[ , dd])^2 ))
#' })
#' z.mgks <- trueF(X)
#' dat <- data.frame(y = z.mgks*cos(z.mgks/2)+rnorm(n))
#' dat$Xks <- data_mgks(y = z, dist = dist)
#' str(dat)
#' head(dat)
#' 
#' ### If the kernel function is isotropic or the distance cannot be written simply 
#' # as sum of the distances elements, a single distance matrices can be passed.
#' dist <- dist[[1]] + dist[[2]]
#' dat <- data.frame(y = z.mgks*cos(z.mgks/2)+rnorm(n))
#' dat$Xks <- data_mgks(y = z, dist = dist)
#' str(dat)
#' head(dat)
#' @export
data_mgks <- function(y, dist){
  pp <- ncol(y)
  if(!is.list(dist)) dist <- list(dist)
  pdist <- sapply(dist, ncol)
  if(any(pdist != pp)) stop("All distance matrices must have the same number of columns as `y`.")
  
  res <- cbind(y, do.call("cbind", dist))
  colnames(res) <- c(rep("y", pp), rep(paste0("d", 1:length(dist)), each = pp))
  res
}









