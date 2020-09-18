########
# Diagonalise penalty
########
.diagPen <- function(X, S, r)
{
  d <- ncol(X)
  
  B <- .getBmatrix(P = S, r = r)
  
  # Reparametrise so that penalty is diagonal
  X <- X %*% B
  
  out <- list("X" = X, 
              "S" = diag(as.numeric(1:d <= r), d, d), 
              "rank" = r)
  
  return( out )
 
} 