#'
#' Multivariate Gaussian kernel smoother and its derivatives
#' 
#' @name mgks
#' @rdname mgks
#' @export mgks
#' @examples
#' ### Example
#' # Simulate values of 2D function
#' n <- 1000
#' X0 <- cbind(runif(n, -1, 1), runif(n, -4, 4))
#' trueF <- function(x) 3 * x[ , 1] + x[ , 2]^2
#' y <- trueF(X0)
#' 
#' # Evaluate kernel smooth on a grid
#' ngr <- 40
#' xseq1 <- seq(-1, 1, length.out = ngr)
#' xseq2 <- seq(-4, 4, length.out = ngr)
#' X <- as.matrix(expand.grid(xseq1, xseq2))
#' 
#' # Here beta[j] is log(1 / sigma[j]) where sigma is standard deviation
#' beta <- c(log(20), log(5))
#' fit <- mgks(y = y, X = X, X0 = X0, beta = beta, deriv = 0)
#' 
#' par(mfrow = c(2, 1))
#' image(x = xseq1, y = xseq2, z = matrix(fit$d0, ngr, ngr), main = "Estimated")
#' image(x = xseq1, y = xseq2, z = matrix(trueF(X), ngr, ngr), main = "Truth")
#' 
#' ### Checking derivatives by finite differences
#' 
#' # Simulate values of 3D function
#' n <- 1000
#' X0 <- cbind(runif(n, -1, 1), runif(n, -4, 4), runif(n, -2, 2))
#' trueF <- function(x) 3 * x[ , 1] + x[ , 2]^2 - 0.5 * x[ , 3]^2
#' y <- trueF(X0)
#' # Define grid over which derivative will be evaluated
#' ngr <- 5
#' xseq1 <- seq(-1, 1, length.out = ngr)
#' xseq2 <- seq(-4, 4, length.out = ngr)
#' xseq3 <- seq(-2, 2, length.out = ngr)
#' X <- as.matrix(expand.grid(xseq1, xseq2, xseq3))
#' 
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(d0 = function(param){
#'   sum(mgks(y = y, X = X, X0 = X0, beta = param, deriv = 0)$d0)
#' },
#' d1 = function(param){
#'   tmp <- colSums(mgks(y = y, X = X, X0 = X0, beta = param, deriv = 1)$d1)
#'   as.list(tmp)
#' },
#' d2 = function(param){
#'   tmp <- colSums(mgks(y = y, X = X, X0 = X0, beta = param, deriv = 2)$d2)
#'   as.list(tmp)
#' },
#' d3 = function(param){
#'   tmp <- colSums(mgks(y = y, X = X, X0 = X0, beta = param, deriv = 3)$d3)
#'   as.list(tmp)
#' })
#' 
#' # Check derivatives up to order 3
#' beta <- c(log(10), log(5), log(5))
#' (tmp <- check_deriv(obj = obj, param = beta, ord = 1:3))
#' 
#' par(mfrow = c(2, 2))
#' for(ii in 1:3){
#'  plot(tmp[[ii]][ , 2], tmp[[ii]][ , 1] - tmp[[ii]][ , 2],
#'       main = paste0("mgks deriv ", ii), ylab = "Error",
#'       xlab = "Deriv value")
#'  abline(h = 0)
#' }
#' 
mgks <- function(y, X, X0, beta, deriv = 0){
  
  if( is.matrix(y) ){ y <- as.vector(y) }
  
  n <- nrow(X)
  n0 <- nrow(X0)
  d <- ncol(X0)
  p <- length(beta)
  stopifnot( length(y) == n0 )
  
  # Under Gaussian log-kernel (x - x0)^2 / sigma^2 = (x - x0)^2 * w hence w = 1/sigma^2 
  w <- exp(2*beta)
  
  tX0 <- t(X0)
  
  d0 <- numeric(n)
  d1 <- d2 <- d3 <- NULL
  if(deriv){ # Allocate storage
    d1 <- matrix(nrow = n, ncol = p)
    if(deriv>1){
      DaDb <- matrix(nrow = n0, ncol = p)
      D2alD2b <- matrix(nrow = n0, ncol = p*(p+1)/2)
      DaDb_DlkDb <- numeric(p*(p+1)/2)
      d2 <- matrix(nrow = n, ncol = p*(p+1)/2)
      if(deriv>2){
        ind2 <- matrix(nrow = p, ncol = p)
        zz <- 1
        for(ir in 1:p){ 
          for (ic in ir:p){ 
            ind2[ir, ic] <- zz
            zz <- zz + 1
          }
        }
        DaDbKK_DlkDb_c <- matrix(nrow = n0, ncol = choose(p + 2, p-1))
        D3alD3b <- matrix(nrow = n0, ncol = choose(p + 2, p-1))
        d3 <- matrix(nrow = n, ncol = choose(p + 2, p-1))
      }
    }
  }
  
  for( ii in 1:n ){
    
    xi <- X[ii, ]
    
    # dist[j, k] = (X0[j, k] - xi[k])^2 * w[k] for j = 1, ..., n0 and k = 1, ..., d
    dist <- - t((tX0 - xi)^2 * w)
    
    # Vector of log-kernels logK[j] = sum_k (X0[j, k] - xi[k])^2 * w[k] for j = 1, ..., n0 
    logK <- rowSums( dist )
    
    # sum exp trick
    mx <- max(logK)
    al <- exp(logK-mx) / sum(exp(logK-mx))
    
    d0[ii] <- sum(al * y)
    
    # Derivatives of g w.r.t. beta
    if( deriv ){
      
      DlkDb <- 2 * dist    # 2 * dist[kk, jj] = D logK[kk] / D beta[jj]  
      
      sum_al_DlkDb <- colSums(al * DlkDb)
      DlkDb_c <- t(t(DlkDb) - sum_al_DlkDb)
      
      # DaDb[i,j] = Dalpha[i] / D beta[j]
      for(jj in 1:p){
        DaDbJJ <- al * DlkDb_c[ , jj]
        d1[ii, jj] <- sum(y * DaDbJJ)
        # Store stuff needed for higher derivs
        if(deriv > 1){ DaDb[ , jj] <- DaDbJJ }
      }
      
      if( deriv > 1){
        
        zz <- 1
        for(jj in 1:p){
          for(kk in jj:p){
            DaDb_DlkDb[zz] <- sum(DaDb[ , kk] * DlkDb[ , jj])
            DaDbKK_DlkDb_cJJ <- DaDb[ , kk] * DlkDb_c[ , jj]
            D2alD2b[ , zz] <- DaDbKK_DlkDb_cJJ - al * DaDb_DlkDb[zz] 
            if(jj == kk){ # Diagonal entries need extra term
              D2alD2b[ , zz] <- D2alD2b[ , zz] + al * 2 * DlkDb_c[ , jj]
            }
            # Store stuff needed for higher derivs
            if( deriv > 2 ){ DaDbKK_DlkDb_c[ , zz] <- DaDbKK_DlkDb_cJJ }
            d2[ii, zz] <- sum(y * D2alD2b[ , zz])
            zz <- zz + 1
          }
        }
        
        if( deriv > 2){
          
          zz <- 1
          for(jj in 1:p){
            for(kk in jj:p){
              for(ll in kk:p){
                D3alD3b[ , zz] <- D2alD2b[ , ind2[kk, ll]] * DlkDb_c[ , jj] - 
                  DaDb[ , kk] * DaDb_DlkDb[ind2[jj, ll]] -
                  DaDb[ , ll] * DaDb_DlkDb[ind2[jj, kk]] - 
                  al * sum(D2alD2b[ , ind2[kk, ll]] * DlkDb[ , jj])
                if(jj == kk){ 
                  D3alD3b[ , zz] <- D3alD3b[ , zz] + 2 * DaDbKK_DlkDb_c[ , ind2[jj, ll]] - 
                    al * 2 * DaDb_DlkDb[ind2[jj, ll]]
                }
                if(jj == ll){
                  D3alD3b[ , zz] <- D3alD3b[ , zz] + 2 * DaDbKK_DlkDb_c[ , ind2[jj, kk]] - 
                    al * 2 * DaDb_DlkDb[ind2[jj, kk]] 
                  if(jj == kk){
                    D3alD3b[ , zz] <- D3alD3b[ , zz] + al * 4 * DlkDb_c[ , jj]
                  }
                }
                d3[ii, zz] <- sum(y * D3alD3b[ , zz])
                zz <- zz + 1
              }
            }
          }
          
        }
      }
      
    }
    
  }
  
  return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
  
}






