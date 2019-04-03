#'
#' Creating a Sinh-Arsinh (shash) family with stack effects
#' 
#' @description This function creates an object representing a SHASH distribution with stack effects
#' @name createSHwithStack
#' @param vector of regression coefficients beta Z.
#' @theta vector of other parameters c(tau, eps phi) of the SHASH distribution
#' @X Model matrix for the convex combination
#' @Z Model matrix
#' @y Response variable
#' @rdname createSHwithStack
#' @export createSHwithStack
#' 
createSHwithStack <- function(param, theta, X, Z, y, S, lam) {
  
  # Loss and derivatives wrt alpha
  derObj <- function(param, deriv = 3) {
    
    # Derivative of the parameter with respect to param, when using identity link
    
    # MODEL: 
    # y_i ~ SHASH(mu_i, tau, eps, phi)
    # tau = 1
    # eps = 1
    # phi = 1
    # h(mu_i) = mu_i = eta_i (identity link)
    # eta_i = X_i' a
    # a = a(nu) (multinomial-like reparametrization)
    # nu_k = Z_i' param_k    k=1,...,K-1
    
    p <- ncol(Z)
    n <- length(y)
    K <- ncol(X)
    param <- matrix(param, nrow = p, ncol = K - 1)
    P <- length(param)
    nu <- Z %*% param
    a <- cbind(1, exp(nu)) / (1 + rowSums(exp(nu)))
    eta <- rowSums(a * X)
    mu <- eta # With identity link mu=eta
    pars <- cbind(mu, theta[1], theta[2], theta[3])
    
    # l wrt mu (actually eta, identity link)
    objSH <- createSH(y = y)$derObj(param = pars, deriv = deriv) 
    
    l <- sum(objSH$d0())
    
    if (deriv > 0) {
      
      objStack <- createStackEffect(X)$derObj(nu, deriv = deriv)
      objEtaNu <- stack.etanu(en = objStack$d1(), 
                              le = objSH$d1(SUM = FALSE)[[1]], 
                              deriv = 1)
      ln <- objEtaNu$ln
      l1 <- as.vector(t(Z) %*% ln)
      
      if (deriv > 1) {
        
        objEtaNu <- stack.etanu(en = objStack$d1(), 
                                le = objSH$d1(SUM = FALSE)[[1]],
                                enn = objStack$d2(), 
                                lee = objSH$d2(SUM = FALSE)[[1]], 
                                deriv = 2)
        lnn <- objEtaNu$lnn
        i2 <- trind.generator(K - 1)$i2
        lbb <- lapply(1:(K - 1), function(x) vector("list", K - 1))
        for (rr in 1:(K - 1)) for (ss in rr:(K - 1)) {
          lbb[[rr]][[ss]] <- t(Z * lnn[, i2[rr, ss]]) %*% Z
          if (ss > rr) lbb[[ss]][[rr]] <- lbb[[rr]][[ss]]
        }
        lbb <- lapply(lbb, function (x) do.call(rbind, x))
        l2 <- matrix(do.call(cbind, lbb), nrow = P)
        
        if (deriv > 2) {
          
          objEtaNu <- stack.etanu(en = objStack$d1(), 
                                  le = objSH$d1(SUM = FALSE)[[1]],
                                  enn = objStack$d2(), 
                                  lee = objSH$d2(SUM = FALSE)[[1]],
                                  ennn = objStack$d3(), 
                                  leee = objSH$d3(SUM = FALSE)[[1]], 
                                  deriv = 3)
          lnnn <- objEtaNu$lnnn
          i3 <- trind.generator(K - 1)$i3
          lbbb <- lapply(1:(K - 1), function(x) {
            lapply(1:(K - 1), function(y) vector("list", K -1))
          })
          for (rr in 1:(K - 1)) for (ss in 1:(K - 1)) for (tt in 1:(K - 1)) {
            lbbb[[rr]][[ss]][[tt]] <- array(0, dim = c(p, p, p))
            for (jj in 1:p) for (kk in 1:p) for (ll in 1:p) {
              lbbb[[rr]][[ss]][[tt]][jj, kk, ll] <- 
                sum(lnnn[, i3[rr, ss, tt]] * Z[, jj] * Z[, kk] * Z[, ll])
            }
          }
          for (rr in 1:(K - 1)) for (ss in 1:(K - 1)) {
            lbbb[[rr]][[ss]] <- do.call(abind, list(lbbb[[rr]][[ss]], along = 3))
          }
          for (rr in 1:(K - 1)) {
            lbbb[[rr]] <- do.call(abind, list(lbbb[[rr]], along = 2))
          }
          lbbb <- do.call(abind, list(lbbb, along = 1))
          
          lbbbMat <- list()
          coun <- 1
          for (jj in 1:P) for (kk in jj:P) for (ll in kk:P) {
            lbbbMat[[coun]] <- lbbb[jj, kk, ll]
            coun <- coun + 1
          }
          l3 <- do.call(c, lbbbMat)
          
        }
      }
    }
    
    # Provide accessors to derivatives
    d0 <- function(){
      return( l )
    }
    d1 <- function(){
      if(deriv < 1) { stop("deriv < 1") }
      return( l1 )
    }
    d2 <- function(){
      if(deriv < 2) { stop("deriv < 2") }
      return( l2 )
    }
    d3 <- function(){
      if(deriv < 3) { stop("deriv < 3") }
      return( l3 )
    }
    
    return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
    
  }
  
  initialize <- function(param, theta, X, Z, y, S, lam) {
    return( createSHwithStack(param, theta, X, Z, y, S, lam) )
  }
  
  return( list("derObj" = derObj, "initialize" = initialize) )
  
}
