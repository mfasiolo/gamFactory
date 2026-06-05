#'
#' Process chicago data from gamair package
#' 
#' @name chicago_data
#' @rdname chicago_data
#' @examples 
#' dat <- chicago_data()
#' dim(dat)
#' str(dat)
#' @export
#' 
#' 

# ============= previous version =============
# chicago_data <- function(si_nexp = FALSE, ...){
#   
#   if( !require(gamair) ){
#     message("Please install the package gamair")
#     return(NULL)
#   }
#   data(chicago, package = "gamair")
#   
#   my_data <- data.frame(death = chicago$death, time = chicago$time)
#   my_data$X <-  as.matrix(chicago[ , c("pm10median", "o3median", "so2median")])
#   my_data$X <- t( t(my_data$X) + apply(my_data$X, 2, function(x) abs(min(x, na.rm = TRUE)))) 
#   my_data$X <- sqrt(1+my_data$X)
#   colnames(my_data$X) <- c("pm10", "o3", "so2")
#   
#   my_data <- na.omit(my_data)
#   
#   # Lagged pollutants for expsmooth
#   my_data$pm_10_lag <- cbind(c(NA, my_data$X[-nrow(my_data), 1]), 
#                              1, c(NA, diff(my_data$time)))
#   colnames(my_data$pm_10_lag) <- c("y", "x", "x")
#   
#   my_data$o3_lag <- cbind(c(NA, my_data$X[-nrow(my_data), 2]), 1, 
#                           c(NA, diff(my_data$time)))
#   colnames(my_data$o3_lag) <- c("y", "x", "x")
#   
#   # Lagged pollutants for mgks
#   n_lag <- 50
#   pm_10_lag_2 <- lapply((n_lag + 1):nrow(my_data), function(ii){
#     kk <- (ii-n_lag):(ii-1)
#     return( c(my_data$X[kk, 1], my_data$time[ii] - my_data$time[kk]) ) 
#   }
#   )
#   my_data$pm_10_lag_2 <- rbind(matrix(NA, n_lag, n_lag*2), do.call("rbind", pm_10_lag_2))
#   colnames(my_data$pm_10_lag_2) <- c(rep("y", n_lag), rep("d1", n_lag))
#   head(my_data$pm_10_lag_2)
#   
#   o3_lag_2 <- lapply((n_lag + 1):nrow(my_data), function(ii){
#     kk <- (ii-n_lag):(ii-1)
#     return( c(my_data$X[kk, 2], my_data$time[ii] - my_data$time[kk]) ) 
#   }
#   )
#   my_data$o3_lag_2 <- rbind(matrix(NA, n_lag, n_lag*2), do.call("rbind", o3_lag_2))
#   colnames(my_data$o3_lag_2) <- c(rep("y", n_lag), rep("d1", n_lag))
#   head(my_data$o3_lag_2)
#   
#   my_data <- na.omit(my_data)
#   
#   if (isTRUE(si_nexp)) {
#     my_data <- pre_si_nexp_chicago(dat = my_data, ...)
#   }
#   
#   return(my_data)
# } 

# ============= new version (without structure fun)=============
chicago_data <- function() {
  if( !require(gamair) ){
    message("Please install the package gamair")
    return(NULL)
  }
  data(chicago, package = "gamair")
  
  
  my_data <- data.frame(death = chicago$death, time = chicago$time)
  my_data$X <-  as.matrix(chicago[ , c("pm10median", "o3median", "so2median")])
  my_data$X <- t( t(my_data$X) + apply(my_data$X, 2, function(x) abs(min(x, na.rm = TRUE)))) 
  my_data$X <- sqrt(1+my_data$X)
  colnames(my_data$X) <- c("pm10", "o3", "so2")
  my_data <- na.omit(my_data)
  
  # Lagged pollutants for expsmooth
  my_data$pm_10_lag <- c(NA, my_data$X[-nrow(my_data), 1])
  my_data$o3_lag <- c(NA, my_data$X[-nrow(my_data), 2])
  my_data$time_gap <- c(NA, diff(my_data$time))
  
  # ===========================================================
  # ❗️Lagged pollutants for mgks(may need to be modified)
  n_lag <- 50
  pm_10_lag_2 <- lapply((n_lag + 1):nrow(my_data), function(ii){
    kk <- (ii-n_lag):(ii-1)
    return( c(my_data$X[kk, 1], my_data$time[ii] - my_data$time[kk]) )
  }
  )
  my_data$pm_10_lag_2 <- rbind(matrix(NA, n_lag, n_lag*2), do.call("rbind", pm_10_lag_2))
  colnames(my_data$pm_10_lag_2) <- c(rep("y", n_lag), rep("d1", n_lag))
  head(my_data$pm_10_lag_2)
  
  o3_lag_2 <- lapply((n_lag + 1):nrow(my_data), function(ii){
    kk <- (ii-n_lag):(ii-1)
    return( c(my_data$X[kk, 2], my_data$time[ii] - my_data$time[kk]) )
  }
  )
  my_data$o3_lag_2 <- rbind(matrix(NA, n_lag, n_lag*2), do.call("rbind", o3_lag_2))
  colnames(my_data$o3_lag_2) <- c(rep("y", n_lag), rep("d1", n_lag))
  head(my_data$o3_lag_2)
  # ===========================================================

  
  # because the first observation was filled by NA
  # not sure if this is necessary
  # ❗️don't remove na: results from fit_si is better (similar to old version)
  # remove na: Results from fit_exp is more like mgks
  my_data <- na.omit(my_data)

  return(my_data)
}
