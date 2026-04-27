#'
#' Compute nested effect's p-value 
#' 
#' @noRd
.sumChi2 <- function(ss, V, Jac) {
  
  Vf <- unname(colSums(t(Jac)* tcrossprod(V, Jac)))
  ss <- ss/Vf
  Jac <- Jac/Vf
  
  # Compute the sum of squared fitted values
  t_y <- sum(ss^2)
  
  # Compute the QR decomposition of the jacobian matrix
  qrx <- qr(Jac, tol = 0)
  R <- qr.R(qrx)
  
  # Compute the covariance matrix V
  V <- R %*% V %*% t(R)
  V <- (V + t(V))/2
  
  # Eigen decomposition of the covariance matrix
  ed <- eigen(V, symmetric = TRUE, only.values = TRUE)
  
  # Compute the p-value using the weighted sum of Chi-square distribution
  pval <- psum.chisq(q = t_y, lb = ed$values, df = rep(1, length(ed$values)))
  
  list(stat=t_y,pval=min(1,pval))
}