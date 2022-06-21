#########################
#' Linearly extrapolates a univariate spline basis
#' 
#' @param x the vector of points at which basis is evaluated.
#' @param b a list where \code{b[[1]]} contains the basis evaluated at x, 
#'          \code{b[[2]]} the matrix of first derivatives (w.r.t. x) and so on.
#' @param th the two boundaries beyond which the linear extrapolation occurs.
#' @param Xbo a 2 x p matrix, where p is the number of basis functions. It contains
#'            the basis function evaluated at \code{th[1]} and \code{th[2]}.
#' @param Xbo1 same as \code{Xbo} but this is the first derivative of the 
#'             basis function.
#' @param method if set to \code{"simple"} higher derivatives are set to 0 beyond 
#'              \code{range(th)}. If set to \code{"smooth"} higher derivatives 
#'              decay to zero at rate r as we move away from 
#'              (\code{th[1]},\code{th[2]}).
#'
#' @return A list where \code{b[[1]]} contains the linearly extrapolated 
#'         basis functions evaluated at \code{x}, \code{b[[2]]} the matrix of 
#'         first derivatives (w.r.t. x) and so on.
#' 
#' @name linextr
#' @rdname linextr
#' @export linextr
#' @examples
#' library(gamFactory)
#' 
#' # Set up basis
#' xseq <- seq(-4, 4, length.out = 1e3) 
#' k <- 20
#' m <- 4
#' ko <- k + m + 2
#' knots <- qunif((1:ko)/(ko+1), -5, 5)
#' 
#' mybasis <- basis_bspline(knots = knots, m = m) 
#' 
#' # Evaluate design matrix and derivatives 
#' X <- basis_bspline(knots = knots, m = m)$evalX(x = xseq, deriv = 3)
#' 
#' # Get limits of "inner" knots and derivatives of basis at the extremal knots
#' krange <- mybasis$krange
#' X_th <- basis_bspline(knots = knots, m = m)$evalX(x = krange, deriv = 4)
#' 
#' # Linearise basis beyond krange
#' X_ext1 <- linextr(x = xseq, b = X, th = krange, 
#'                   Xbo = X_th$X0, Xbo1 = X_th$X1, method = "simple")
#' 
#' # Smoothly linearise basis beyond krange
#' X_ext2 <- linextr(x = xseq, b = X, th = krange, 
#'                   Xbo = X_th$X0, Xbo1 = X_th$X1, method = "smooth")
#' 
#' # Plot effect with random coefficients
#' betas <- runif(ncol(X$X0))
#' par(mfrow = c(1, 1))
#' plot(xseq, X_ext1$X0%*%betas, type = 'l', col = 1) # Linearly extrapolated
#' lines(xseq, X_ext2$X0%*%betas, type = 'l', col = 4)# Smoothly linearly extrap
#' lines(xseq, X$X0%*%betas, col = 2)                 # Original
#' 
#' # Above we see some small bumps just beyond krange (look at discrepancy
#' # between blue and red, you might to run this a few times to see it).
#' # We can reduce this by setting to zero a few derivatives of order > 1 
#' # at krange. For i = 1 and 2, below we impose that
#' # f''(krange[i]) = f'''(krange[i]) = f''''(krange[i]) =  0  
#' 
#' # We first constrain the original basis by reparametrising it
#' library(MASS)
#' NS <- Null( t(rbind(X_th$X2, X_th$X3, X_th$X4)) )
#' X_con <- X
#' X_con$X0 <- X$X0 %*% NS
#' X_con$X1 <- X$X1 %*% NS
#' X_con$X2 <- X$X2 %*% NS
#' X_con$X3 <- X$X3 %*% NS
#' 
#' # Linearise reparametrisd basis beyond krange
#' X_ext1 <- linextr(x = xseq, b = X_con, th = krange, 
#'                   Xbo = X_th$X0%*%NS, Xbo1 = X_th$X1%*%NS, method = "simple")
#' 
#' # Smoothly linearise basis beyond krange
#' X_ext2 <- linextr(x = xseq, b = X_con, th = krange, 
#'                   Xbo = X_th$X0%*%NS, Xbo1 = X_th$X1%*%NS, method = "smooth")
#' 
#' # Plot effect with random coefficients
#' betas <- runif(ncol(X_con$X0))
#' par(mfrow = c(1, 1))
#' plot(xseq, X_ext1$X0%*%betas, type = 'l', col = 1) #Linearly extrapolated
#' lines(xseq, X_ext2$X0%*%betas, type = 'l', col = 4)#Smoothly linearly extrap
#' lines(xseq, X_con$X0%*%betas, col = 2)             #Original (reparametrised)
#' # No appreciable difference between standard (black) and smooth (blue)
#' # extrapolation
#' 
#' # Checking whether the derivatives of the extrapolated basis are right
#' wrap_basis <- function(x){
#'   X <- basis_bspline(knots = knots, m = m)$evalX(x = x, deriv = 3)
#'   X$X0 <- X$X0 %*% NS
#'   X$X1 <- X$X1 %*% NS
#'   X$X2 <- X$X2 %*% NS
#'   X$X3 <- X$X3 %*% NS
#'   X <- linextr(x = x, b = X, th = krange, 
#'                Xbo = X_th$X0%*%NS, Xbo1 = X_th$X1%*%NS, method = "smooth")
#'   return(X) 
#' }
#' 
#' l0 <- wrap_basis(x = xseq)
#' lm <- wrap_basis(x = xseq - 1e-5)
#' lp <- wrap_basis(x = xseq + 1e-5)
#' 
#' # It seems so!
#' par(mfrow = c(2, 2))
#' plot(l0$X1, (lp$X0 - lm$X0) / 2e-5)
#' abline(0, 1, col = 2)
#' plot(l0$X2, (lp$X1 - lm$X1) / 2e-5)
#' abline(0, 1, col = 2)
#' plot(l0$X3, (lp$X2 - lm$X2) / 2e-5)
#' abline(0, 1, col = 2)
#' 
#' # Note that with smooth extrapolation, derivative up to order 3 are continuous
#' # (assuming that the derivative of the original basis are continuous!)
#' xseq <- seq(-4, 4, length.out = 1e3) 
#' l0 <- wrap_basis(x = xseq)
#' 
#' par(mfrow = c(2, 2))
#' matplot(xseq, l0$X0, type = 'l') # 1st derivative w.r.t. x
#' matplot(xseq, l0$X1, type = 'l') # 1st derivative w.r.t. x
#' matplot(xseq, l0$X2, type = 'l') # 2nd derivative w.r.t. x
#' matplot(xseq, l0$X3, type = 'l') # 3rd derivative w.r.t. x
#' 
#' # Zooming in where the linearisation occurs
#' xseq <- seq(krange[2]-0.3, krange[2]+0.3, length.out = 1000)
#' l0 <- wrap_basis(x = xseq)
#' 
#' par(mfrow = c(2, 2))
#' matplot(xseq, l0$X0, type = 'l') # 1st derivative w.r.t. x
#' matplot(xseq, l0$X1, type = 'l') # 1st derivative w.r.t. x
#' matplot(xseq, l0$X2, type = 'l') # 2nd derivative w.r.t. x
#' matplot(xseq, l0$X3, type = 'l') # 3rd derivative w.r.t. x
#' abline(v = krange)
#'
linextr <- function(x, b, th, Xbo, Xbo1, method = "smooth", r = NULL){
  
  deriv <- length(b) - 1
  
  # Find x's falling outside (th[1], th[2])
  d <- ncol(b$X0)
  kl <- which(x < th[1])
  ku <- which(x > th[2])
  kk <- c(kl, ku)
  
  # If any x outside, get linear approximation to the basis
  if( length(kk) > 0 ) {
    # Create slopes depending on whether we are beyond lower or upper boundary [a, b]
    del <- c()
    Xth <- Xth1 <- matrix(NA, ncol = ncol(Xbo), nrow = 0)
    if( length(kl) > 0 ){ 
      del <- c(del, x[kl] - th[1])
      Xth <- rbind(Xth, matrix(Xbo[1, ], length(kl), d, byrow = TRUE))
      Xth1 <- rbind(Xth1, matrix(Xbo1[1, ], length(kl), d, byrow = TRUE))
    }
    if( length(ku) > 0 ){ 
      del <- c(del, x[ku] - th[2])
      Xth <- rbind(Xth, matrix(Xbo[2, ], length(ku), d, byrow = TRUE))
      Xth1 <- rbind(Xth1, matrix(Xbo1[2, ], length(ku), d, byrow = TRUE))
    }
    D
    # Simple method: compute value of 1st order Taylor expansion evaluatated at each point.
    # Further derivatives are set to zero.
    if(method == "simple"){
      b$X0[kk, ] <- Xth + Xth1 * del
      if(deriv){
        b$X1[kk, ] <- Xth1
        if(deriv > 1){
          for(id in 3:(deriv+1))
            b[[id]][kk, ] <- 0
        }
      }
    }
    
    # Smooth decay method: exponential decay of higher derivatives beyond thresholds [a, b].
    # Decay is exp(-r*|x-c|^3) where c = {a or b}. 
    if(method == "smooth"){
      
      # Set r = K / |b - a|^3 for some K > 0 to get decay exp(-K * (|x-c|/(b-a))^3 )
      if( is.null(r) ) { r <- 1e4 / diff(range(th))^3 } 
      
      deriv <- length(b) - 1
      if(deriv > 3){ stop("Cannot handle derivatives of order higher than 3.") }
      
      a <- exp(-r*abs(del)^3)
      XL <- Xth + Xth1 * del
      X0 <- b$X0[kk, ]
      b$X0[kk, ] <- a * X0 + (1-a) * XL
      
      if(deriv){
        r3del2 <- sign(del)*3*r*del^2
        X1 <- b$X1[kk, ]
        X0mL <- X0-XL
        A <- (X1 - r3del2*X0mL)
        b$X1[kk, ] <- a*A + (1-a)*Xth1
        if(deriv > 1){
          X2 <- b$X2[kk, ]
          X1mL1 <- X1-Xth1
          L1mA <- Xth1 - A
          A1 <- X2 - r3del2/del*(2*X0mL + del*X1mL1)
          b$X2[kk, ] <- a*(r3del2*L1mA + A1)
          if(deriv > 2){
            X3 <- b$X3[kk, ]
            A2 <- X3 - sign(del)*3*r*(2*X0mL+del*(4*X1mL1+del*X2))
            b$X3[kk, ] <- -r3del2*b$X2[kk, ] + a*(r3del2/del*(2*L1mA-del*A1)+A2)
          }
        }
      }
    }
  }
  
  
  return(b)
  
}