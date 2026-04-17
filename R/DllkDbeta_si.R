#'
#' Derivatives of single index effects
#' 
#' @rdname DllkDbeta.si
#' @export DllkDbeta.si
#' @export
#'
DllkDbeta.si <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv == 0 ){ return( list() ) }
  
  if( is.null(param) ){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  }
  
  # Need to update the object
  if( is.null(o$param) || !identical(param, o$param) || o$deriv < deriv){
    o <- o$eval(param = param, deriv = deriv)
  }
  
  # Extract stuff from object for convenience
  alpha <- param[ 1:o$na ]
  beta <- param[ -(1:o$na) ]
  
  na <- length( alpha )
  nb <- length( beta ) 
  
  Xi <- o$store$Xi
  
  # [新增] 提取正约束开关和 a0。注意这里的 o 结构中 si 信息通常在 o$xt$si
  positive_si <- isTRUE(o$xt$si$positive_si)
  if(positive_si) {
    a0 <- o$xt$si$a0
    if(is.null(a0)) a0 <- rep(0, na)
    
    # 计算指数映射因子
    exp_alpha <- exp(alpha + a0)
    
    # 核心步骤：直接通过缩放 Xi 来吸收一阶链式法则的 exp() 乘子
    Xi <- t(t(Xi) * exp_alpha)
  }
  
  # Outer design matrix and its derivatives
  X <- o$store$X0; X1 <- o$store$X1; X2 <- o$store$X2; X3 <- o$store$X3
  
  f <- o$f
  
  der1 <- der2 <- der3 <- der4 <- NULL
  
  f1 <- drop( X1 %*% beta )
  le <- llk$d1
  lg <- le * f1
  
  ll_a <- t(Xi) %*% lg
  ll_b <- t(X) %*% le
  
  der1 <- c(ll_a, ll_b)
  
  if( deriv >  1){
    f2 <- drop( X2 %*% beta )
    lee <- llk$d2
    lgg <- le * f2 + lee * f1^2
    leg <- lee * f1
    
    ll_aa <- t(Xi) %*% (lgg * Xi) 
    
    # [新增] 二阶导数的海森矩阵修正
    if (positive_si) {
      # 根据链式法则：d^2(l)/d(alpha_inner)^2 必须加上一阶导数 lg * d^2(xa)/d(alpha_inner)^2
      # 因为 d^2(exp)/dx^2 还是 exp，所以这个修正项完美等同于对角线化的一阶导数 ll_a
      ll_aa <- ll_aa + diag(as.vector(ll_a))
    }
    
    ll_bb <- t(X) %*% (lee * X)    # Same as t(X) %*% diag(le2) %*% X
    ll_ba <- t(X) %*% (leg * Xi) + t(X1) %*% (le * Xi)
    
    der2 <- rbind(cbind(ll_aa, t(ll_ba)), 
                  cbind(ll_ba, ll_bb))
    
    der3 <- list()
    if( deriv >  2  ){
      f3 <- drop( X3%*%beta )
      leee <- llk$d3
      lggg <- le * f3 + 3 * lee * f1 * f2 + leee * f1^3
      leeg <- leee * f1
      legg <- leee * f1^2 + lee * f2
      
      coun <- 1
      for(jj in 1:na){      # AXX
        XJ <- Xi[ , jj]
        for(kk in jj:na){   # AAX
          XJK <- XJ * Xi[ , kk]
          for(ll in kk:na){ # AAA
            der3[[coun]] <- sum(lggg * XJK * Xi[ , ll])
            coun <- coun + 1
          }
          for(ll in 1:nb){  # AAB
            der3[[coun]] <- sum((legg*X[,ll]+2*leg*X1[,ll]+le*X2[,ll]) * XJK)
            coun <- coun + 1
          }
        }
        for(kk in 1:nb){    # ABB
          XJK <- XJ * X[ , kk]
          for(ll in kk:nb){
            der3[[coun]] <- sum(XJ * (leeg*X[,ll]*X[,kk] + lee*(X1[,ll]*X[,kk]+X[,ll]*X1[,kk]) ) )
            coun <- coun + 1
          }
        }
      }
      for(jj in 1:nb){      # BBB
        XJ <- leee * X[ , jj]
        for(kk in jj:nb){
          XJK <- XJ * X[ , kk]
          for(ll in kk:nb){
            der3[[coun]] <- sum( XJK * X[ , ll] ) # Same as t(leee) %*% (X[ , jj] * X[ , kk] * X[ , ll])
            coun <- coun + 1
          }
        }
      }
      der3 <- do.call("c", der3)
      
    } # deriv = 3
  }   #         2
  
  out <- list("d1" = der1, "d2" = der2, "d3" = der3)
  
  return( out )
  
}
