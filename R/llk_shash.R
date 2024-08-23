#'
#' Log-likelihood of the Sinh-Arsinh (shash) distribution
#' 
#' @description Log-likelihood of the Sinh-Arsinh (shash) distribution 
#' @param y a vector of observations.
#' @param phiPen positive coefficient of ridge penalty on phi parameter.
#' @examples 
#' 
#' library(gamFactory)
#' 
#' ##### Simulate data and create object
#' n <- 10
#' param <- c(1, exp(2), 0.5, exp(1))
#' y <- param[1] + (param[4] * param[2]) * sinh((1/param[4]) * asinh(qnorm(runif(n))) + (param[3]/param[4]))
#' 
#' llk_shash(y = y, param = param, deriv = 3)
#' 
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_shash(y = y, param = param, deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_shash(y = y, param = param, deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_shash(y = y, param = param, deriv = 2)$d2))
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_shash(y = y, param = param, deriv = 3)$d3))
#'   })
#' 
#' check_deriv(obj = obj, param = param, ord = 1:3)
#' @export
#'
llk_shash <- function(y, param, deriv = 0, phiPen = 1e-3){
  
  sech <- function(.x){ 1 / cosh(.x) }
  
  .log1pexp <- function(x) 
  {
    indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), right = TRUE, 
                     include.lowest = TRUE)
    kk <- which(indx == 1)
    if (length(kk)) {
      x[kk] <- exp(x[kk])
    }
    kk <- which(indx == 2)
    if (length(kk)) {
      x[kk] <- log1p(exp(x[kk]))
    }
    kk <- which(indx == 3)
    if (length(kk)) {
      x[kk] <- x[kk] + exp(-x[kk])
    }
    return(x)
  }
  
  # Compute sqrt(x^2 + m) when |x| >> 0 and m is reasonably small (e.g. + 1 or - 1)
  .sqrtX2pm <- function(x, m){ 
    x <- abs(x)
    kk <- which( x < 1e8 )
    if( length(kk) ){
      x[kk] <- sqrt(x[kk]^2 + m)
    }
    return(x)
  }
  
  # Compute (a*x^2 + m1) / (x^2 + m2)^2 when |x| >> 0 and m1, m2 are reasonably small (e.g. + 1 or - 1)
  .ax2m1DivX2m2SQ <- function(x, m1, m2, a = 1){
    if(a < 0){ stop("'a' has to be positive")  }
    x <- abs(x)
    kk <- (a * x^2 + m1) < 0
    o <- x * 0
    if( any(kk) ){
      o[kk] <- (a * x[kk]^2 + m1) / (x[kk]^2 + m2)^2
    }
    if( sum(kk) < length(x) ){
      o[!kk] <- ((.sqrtX2pm(sqrt(a)*x[!kk], m1) / .sqrtX2pm(x[!kk], m2)) / .sqrtX2pm(x[!kk], m2))^2
    }
    return(o)
  }
  
  if (is.list(param) ) param <- do.call("cbind", param)
  if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
  if (ncol(param) != 4) stop("Wrong number of parameters provided")
  
  mu  <- param[ , 1, drop = TRUE]
  tau <- param[ , 2, drop = TRUE]
  eps <- param[ , 3, drop = TRUE]
  phi <- param[ , 4, drop = TRUE]  
  
  sig <- exp( tau )
  del <- exp( phi )
  
  # 1) Calculate derivative of likelihood to appropriate order
  z <- (y - mu) / (sig*del)
  
  dTasMe <- del*asinh(z) - eps
  g <- -dTasMe
  CC <- cosh( dTasMe )
  SS <- sinh( dTasMe )
  
  l <- -tau - 0.5*log(2*pi) + log(CC) - 0.5*.log1pexp(2*log(abs(z))) - 0.5*SS^2 - phiPen*phi^2
  
  out <- list()
  out$d0 <- l
  
  if( deriv > 0 )
  {
    zsd <- z*sig*del
    sSp1 <- .sqrtX2pm(z, 1) # sqrt(z^2+1)
    asinhZ <- asinh(z)
    
    ## First derivatives 
    De <- tanh(g) - 0.5*sinh(2*g)
    Dm <- 1/(del*sig*sSp1)*(del*(De)+z/sSp1)
    Dt <- zsd*Dm - 1
    Dp <- Dt + 1 - del*asinhZ*De - 2*phiPen*phi
    out[["d1"]] <- list(Dm, Dt, De, Dp) 
    
    if( deriv > 1 ){
      Dme <- (sech(g)^2 - cosh(2*g)) / (sig*sSp1)
      Dte <- zsd*Dme
      Dmm <- Dme/(sig*sSp1) + z*De/(sig^2*del*sSp1^3) + .ax2m1DivX2m2SQ(z, -1, 1)/(del*sig*del*sig)
      Dmt <- zsd*Dmm - Dm
      Dee <- -2*cosh(g)^2 + sech(g)^2 + 1 
      Dtt <-  zsd*Dmt
      Dep <- Dte - del*asinhZ*Dee
      Dmp <- Dmt + De/(sig*sSp1) - del*asinhZ*Dme
      Dtp <- zsd*Dmp
      Dpp <- Dtp - del*asinhZ*Dep + del*(z/sSp1-asinhZ)*De - 2*phiPen
      out[["d2"]] <- list(l_mm = Dmm, l_mt = Dmt, l_me = Dme, l_mp = Dmp, l_tt = Dtt, 
                          l_te = Dte, l_tp = Dtp, l_ee = Dee, l_ep = Dep, l_pp = Dpp)
      
      if( deriv > 2 ){
        Deee <-  -2*(sinh(2*g)+sech(g)^2*tanh(g))
        Dmee <- Deee/(sig*sSp1)
        Dmme <- Dmee/(sig*sSp1) + z*Dee/(sig*sig*del*sSp1^3)
        Dmmm <- 2*z*Dme/(sig*sig*del*sSp1^3) + Dmme/(sig*sSp1) + 
          .ax2m1DivX2m2SQ(z, -1, 1, 2)*De/(sig^3*del^2*sSp1) + 
          2*(z/sSp1)*.ax2m1DivX2m2SQ(z, -3, 1)/((sig*del)^3*sSp1)
        Dmmt <- zsd*Dmmm - 2*Dmm
        Dtee <- zsd*Dmee
        Dmte <- zsd*Dmme - Dme
        Dtte <- zsd*Dmte
        Dmtt <- zsd*Dmmt - Dmt
        Dttt <- zsd*Dmtt
        Dmep <- Dmte + Dee/(sig*sSp1) - del*asinhZ*Dmee
        Dtep <- zsd*Dmep
        Deep <- Dtee - del*asinhZ*Deee
        Depp <- Dtep - del*asinhZ*Deep + del*( z/sSp1-asinhZ )*Dee
        Dmmp <- Dmmt + 2*Dme/(sig*sSp1) + z*De/(del*sig*sig*sSp1^3) - del*asinhZ*Dmme
        Dmtp <- zsd*Dmmp - Dmp
        Dttp <- zsd*Dmtp
        Dmpp <- Dmtp + Dep/(sig*sSp1) + z^2*De/(sig*sSp1^3) - 
          del*asinhZ*Dmep + del*Dme*(z/sSp1 - asinhZ)
        Dtpp <- zsd*Dmpp
        Dppp <- Dtpp - del*asinhZ*Depp + del*(z/sSp1-asinhZ)*(2*Dep + De) + del*(z/sSp1)^3 * De
        out[["d3"]] <- list(l_mmm = Dmmm, l_mmt = Dmmt, l_mme = Dmme, l_mmp = Dmmp,
                            l_mtt = Dmtt, l_mte = Dmte, l_mtp = Dmtp, l_mee = Dmee,
                            l_mep = Dmep, l_mpp = Dmpp, l_ttt = Dttt, l_tte = Dtte,
                            l_ttp = Dttp, l_tee = Dtee, l_tep = Dtep, l_tpp = Dtpp,
                            l_eee = Deee, l_eep = Deep, l_epp = Depp, l_ppp = Dppp) 
      }
    }
  }
  
  return( out )
  
}