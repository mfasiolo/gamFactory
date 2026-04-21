#'
#' Predict using single index effects
#' 
#' @noRd
.predict.matrix.si <- function(object, data, get.xa = FALSE){
  
  # Need to compute single index vector by projecting inner model matrix on alpha
  si <- object$xt$si
  
  alpha <- si$alpha
  a0 <- si$a0
  
  # Need to subtract colMeans of original data "xm" and rescale using B
  Xi <- t(t(data[[object$term]]) - si$xm)  %*% si$B
  
  positive_si <- isTRUE(si$positive_si) # 安全获取约束开关
  # 根据约束开关计算 xa (预测时的单指数向量)
  if (positive_si) {
    exp_alpha_a0 <- exp(alpha + a0)
    xa <- Xi %*% exp_alpha_a0
  } else {
    xa <- Xi %*% (alpha + a0)
  }
  
  if(get.xa){ 
    # 计算 xa 关于内层参数 alpha 的导数矩阵
    if (positive_si) {
      # 链式法则：d(Xi * exp(alpha + a0)) / d_alpha = Xi * diag(exp(alpha + a0))
      # R 语言里最高效的按列相乘方式是用 t(t(matrix) * vector) 或 sweep
      xa_da <- t(t(Xi) * as.vector(exp_alpha_a0))
    } else {
      xa_da <- Xi
    }
    return(list(xa = xa, xa_da = xa_da))
  }
  
  # Compute outer model matrix
  X0 <- object$xt$basis$evalX(x = xa, deriv = 0)$X0
  
  # Total model matrix is X0 preceded my matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (this is a trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xa
  
  return(Xtot)
}
