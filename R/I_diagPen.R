########
# Diagonalise a penalty
########
# INPUT
# X is the model matrix
# S is the penalty matrix on the coefficients
# r is rank of penalty matrix S
#
# OUTPUT
# Reparametrised versions of B and S, and reparametrisation matrix B
#
.diagPen <- function(X, S, r)
{
  d <- ncol(X)
  
  B <- .getBmatrix(P = S, r = r)
  
  # Reparametrise so that penalty is diagonal
  X <- X %*% B
  
  out <- list("X" = X, 
              "S" = diag(as.numeric(1:d <= r), d, d), 
              "B" = B,
              "rank" = r)
  
  return( out )
 
} 