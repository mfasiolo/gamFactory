##################
#'
#' Fake log-likelihood useful for testing
#'
#' @description A synthetic function (not a real log-likelihood) whose derivatives can be
#'              derived by hand, used to test derivative machinery such as
#'              [`gamFactory::check_deriv`]. Returns
#'              \code{d0(param) = prod_i f_i(param[, i])}, where \code{f_i = cos} if
#'              \code{i} is odd and \code{f_i = sin} if \code{i} is even, together with
#'              its derivatives up to order 3, in the list format used by
#'              [`gamFactory::llk_gaussian`] and friends.
#' @param param a matrix (or list) with one column (element) per dimension.
#' @param deriv integer between 0 and 3 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0}, while 1-3 additionally return
#'              \code{d1}-\code{d3}.
#' @rdname llk_test
#' @export llk_test
#' @examples
#' library(gamFactory)
#' p <- 3
#' param <- matrix(rnorm(p), nrow = 1)
#' llk_test(param = param, deriv = 3)
#'
llk_test <- function(param, deriv){
  
  if( is.vector(param) ) { param <- matrix(param, nrow = 1) }
  
  if( is.list(param) ) { param <- do.call("cbind", param) }
   
  p <- ncol(param)
  
  out <- list()
  
  cosf <- lapply(1:p, function(ii) if(ii%%2) { cos(param[ , ii]) } else { sin(param[ , ii]) }  )
  cosf1 <- lapply(1:p, function(ii) if(ii%%2) { -sin(param[ , ii]) } else { cos(param[ , ii]) }  )
  cosf2 <- lapply(1:p, function(ii) if(ii%%2) { -cos(param[ , ii]) } else { -sin(param[ , ii]) }  )
  cosf3 <- lapply(1:p, function(ii) if(ii%%2) { sin(param[ , ii]) } else { -cos(param[ , ii]) }  )
  
  out$d0 <- Reduce("*", cosf)
  
  if(deriv >= 1){
    out$d1 <- lapply(1:p, function(ii) out$d0 / cosf[[ii]] * cosf1[[ii]])
    
    if( deriv >= 2){
      out$d2 <- list()
      jj <- 1
      for(ir in 1:p){
        for(ic in ir:p){
          if( ir != ic ){
            out$d2[[jj]] <- out$d0 / cosf[[ir]] /  cosf[[ic]] * cosf1[[ir]] * cosf1[[ic]]     
          } else {
            out$d2[[jj]] <- out$d0 / cosf[[ir]] *  cosf2[[ir]] 
          }
          jj <- jj + 1
        }
      }
      
      out$d3 <- list()
      jj <- 1
      for(i1 in 1:p){
        for(i2 in i1:p){
          for(i3 in i2:p){
            if( i1 == i2 && i2 == i3 ){
              out$d3[[jj]] <- out$d0 / cosf[[i1]] * cosf3[[i1]]    
            } else {
              if( i1 == i2 ){
                out$d3[[jj]] <- out$d0 / cosf[[i1]] / cosf[[i3]] * cosf2[[i1]] * cosf1[[i3]]    
              } else {
                if( i2 == i3 ){
                  out$d3[[jj]] <- out$d0 / cosf[[i1]] / cosf[[i2]] * cosf2[[i2]] * cosf1[[i1]]    
                } else {
                  out$d3[[jj]] <- out$d0 / cosf[[i1]] / cosf[[i2]] / cosf[[i3]] * cosf1[[i1]] * cosf1[[i2]] * cosf1[[i3]]  
                }
              }
            }
            jj <- jj + 1    
          }
        }
      }
    }
  }

  return( out )
}


