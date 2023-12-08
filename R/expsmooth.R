#'
#' Exponential smooth and its derivatives
#' 
#' @name expsmooth
#' @rdname expsmooth
#' @export expsmooth
#' @importFrom Rcpp evalCpp
#' @examples
#' n <- 1000
#' xseq <- seq(-4, 4, length.out = n)
#' y <- sin(xseq^2) + rnorm(n, 0, 0.5)
#' plot(xseq, y, col = "grey")
#' 
#' Xi <- cbind(1, xseq, xseq^2)
#' beta <- c(2, 0, -0.2)
#' p <- length(beta)
#' xsm <- expsmooth(y = y, Xi = Xi, beta = beta, deriv = 3)
#' 
#' lines(xseq, xsm$d0, col = 2, lwd = 2)
#' 
#' obj <- list(d0 = function(param){ 
#'   sum(expsmooth(y = y, Xi = Xi, beta = param, deriv = 0)$d0)
#' }, 
#' d1 = function(param){
#'   tmp <- colSums(expsmooth(y = y, Xi = Xi, beta = param, deriv = 1)$d1)
#'   as.list(tmp)
#' }, 
#' d2 = function(param){
#'   tmp <- colSums(expsmooth(y = y, Xi = Xi, beta = param, deriv = 2)$d2)
#'   as.list(tmp)
#' }, 
#' d3 = function(param){
#'   tmp <- colSums(expsmooth(y = y, Xi = Xi, beta = param, deriv = 3)$d3)
#'   as.list(tmp)
#' })
#' 
#' check_deriv(obj = obj, param = beta, ord = 1:3)
#' 
expsmooth <- function(y, Xi, beta, x0 = NULL, deriv = 0){
  
  if( is.matrix(y) ){ y <- as.vector(y) }
  
  if( is.null(x0) ){ x0 <- y[1] }
  
  return( .expsmooth_cpp(y = y, Xi = Xi, beta = beta, x0 = x0, deriv = deriv) )

}

######## OLD R Version
# expsmooth <- function(y, Xi, beta, x0 = NULL, deriv = 0){
#   
#   if( is.matrix(y) ){ y <- as.vector(y) }
#   
#   if( is.null(x0) ){ x0 <- y[1] }
#   
# n <- length(y)
# p <- length(beta)
# 
# w <- plogis( drop(Xi %*% beta) )
# 
# # Perform smoothing
# g <- numeric(n)
# g[1] <- w[1] * x0 + (1-w[1]) * y[1]
# for(ii in 2:n){
#   g[ii] <- w[ii] * g[ii-1] + (1-w[ii]) * y[ii]
# }
# 
# # Derivatives of g w.r.t. beta
# d1 <- d2 <- d3 <- NULL
# if( deriv ){
#   a1 <- w * (1 - w)
#   w1 <- a1 * Xi
#   gmx <- c(x0, g[1:(n-1)]) - y
#   
#   d1 <- matrix(nrow = n, ncol = p)
#   d1[1, ] <- w1[1, ] * gmx[1]
#   for(ii in 2:n){
#     d1[ii, ] <- w1[ii, ] * gmx[ii] + w[ii] * d1[ii-1, ]
#   }
#   
#   if( deriv > 1){
#     a2 <- a1 - 2 * a1 * w
#     w2 <- matrix(nrow = n, ncol = p*(p+1)/2)
#     d2 <- matrix(nrow = n, ncol = p*(p+1)/2)
#     zz <- 1
#     for(jj in 1:p){
#       for(kk in jj:p){
#         w2[ , zz] <- a2 * Xi[,jj] * Xi[,kk]
#         d2[1, zz] <- w2[1, zz] * gmx[1]
#         for(ii in 2:n){
#           d2[ii, zz] <- w2[ii,zz]*gmx[ii] + w1[ii,jj]*d1[ii-1,kk] + 
#             w1[ii,kk]*d1[ii-1,jj] + w[ii]* d2[ii-1, zz]
#         }
#         zz <- zz + 1
#       }
#     }
#     
#     if( deriv > 2){
#       a3 <- a2 - 2*a2*w - 2*a1*a1
#       w3 <- matrix(nrow = n, ncol = choose(p + 2, p-1))
#       d3 <- matrix(nrow = n, ncol = choose(p + 2, p-1))
#       ind <- matrix(nrow = p, ncol = p)
#       zz <- 1
#       for(ir in 1:p){ 
#         for (ic in ir:p){ 
#           ind[ir, ic] <- zz
#           zz <- zz + 1
#         }
#       }
#       zz <- 1
#       for(jj in 1:p){
#         for(kk in jj:p){
#           for(ll in kk:p){
#             w3[ , zz] <- a3 * Xi[,jj] * Xi[,kk] * Xi[,ll] 
#             d3[1, zz] <- w3[1, zz] * gmx[1]
#             for(ii in 2:n){
#               d3[ii,zz] <- w3[ii,zz]*gmx[ii] + 
#                 w2[ii,ind[jj,kk]]*d1[ii-1,ll] + 
#                 w2[ii,ind[jj,ll]]*d1[ii-1,kk] + 
#                 w2[ii,ind[kk,ll]]*d1[ii-1,jj] + 
#                 w1[ii,jj]*d2[ii-1,ind[kk,ll]] + 
#                 w1[ii,kk]*d2[ii-1,ind[jj,ll]] +
#                 w1[ii,ll]*d2[ii-1,ind[jj,kk]] +
#                 w[ii]* d3[ii-1, zz]
#             }
#             zz <- zz + 1
#           }
#         }
#       }
#     }
#   }
# }
# 
# return( list("d0" = g, "d1" = d1, "d2" = d2, "d3" = d3) )





 