#' Process Chicago Air Pollution and Mortality Data
#' 
#' @description Loads, transforms, and structures the Chicago Air Pollution and Mortality Data
#'  from the \code{gamair} package for testing nested effects (\code{\link{s_nest}}).
#' 
#' @details The function executes the following preprocessing steps on the raw data:
#' \enumerate{
#'   \item Pollutant observations (\code{pm10median}, \code{o3median}, \code{so2median})
#'         are shifted by adding the absolute value of their respective minimums to guarantee non-negative values.
#'   \item A square-root transformation is applied via \eqn{\sqrt{1 + X}}.
#'   \item Creates lag 1 vectors for \eqn{PM_{10}} and \eqn{O_3} (\code{pm_10_lag}, \code{o3_lag}) along with a 
#'   vector of time gaps \code{time_gap}, which indicate the time passed between the current observation and the previous.
#'   \item Generates wide matrices containing, for each row, the 
#'         previous \code{n_lag = 50} historical days of pollutant exposure and 
#'         corresponding time distances for localized smoothing operations.
#' }
#' 
#' Final rows containing any \code{NA} values (including the initial 50 rows used 
#' to build history windows and the final row from shifting) are automatically dropped via 
#' \code{na.omit}.
#' 
#' @return A processed \code{data.frame} containing:
#' \itemize{
#'   \item{\code{death}}{ Numeric vector of daily mortality counts.}
#'   \item{\code{time}}{ Numeric vector of sequential time step indices.}
#'   \item{\code{X}}{ A matrix with columns \code{pm10}, \code{o3}, and \code{so2} representing transformed concentrations.}
#'   \item{\code{pm_10_lag / o3_lag}}{ Lag 1 pollutant vectors.}
#'   \item{\code{time_gap}}{ Numeric vector tracking the difference between subsequent time steps.}
#'   \item{\code{pm_10_lag2 / o3_lag2}}{ Matrices of dimension \code{N} by 50 containing rolling historical pollutant entries.}
#'   \item{\code{time_dist}}{ Matrix of dimension \code{N} by 50 tracking historical day gaps relative to the current evaluation row.}
#' }
#' 
#' @name chicago_data
#' @rdname chicago_data
#' @examples 
#' dat <- chicago_data()
#' dim(dat)
#' str(dat)
#' @export
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
  my_data$pm_10_lag <- c(my_data$X[-1, 1], NA)
  my_data$o3_lag <- c(my_data$X[-1, 2], NA)
  my_data$time_gap <- c(NA, diff(my_data$time))
  
  # ===========================================================
  # Lagged pollutants for mgks
  n_lag <- 50
  pm10_yd <- lapply((n_lag + 1):nrow(my_data), function(ii){
    kk <- (ii-n_lag):(ii-1)
    return(list(y = my_data$X[kk, "pm10"], d = my_data$time[ii] - my_data$time[kk]))
  })
  my_data$pm_10_lag2 <- rbind(matrix(NA, n_lag,n_lag), t(sapply(pm10_yd, "[[", "y")))
  colnames(my_data$pm_10_lag2) <- rep("y", n_lag)
  my_data$time_dist <- rbind(matrix(NA, n_lag,n_lag), t(sapply(pm10_yd, "[[", "d")))
  colnames(my_data$pm_10_lag2) <- rep("d", n_lag)
  
  o3_yd <- lapply((n_lag + 1):nrow(my_data), \(ii){ my_data$X[(ii-n_lag):(ii-1), "o3"]})
  my_data$o3_lag2 <- rbind(matrix(NA, n_lag,n_lag), do.call("rbind", o3_yd))
  # ===========================================================
  
  my_data <- na.omit(my_data)
  
  return(my_data)
}