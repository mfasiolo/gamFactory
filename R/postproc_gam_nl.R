#' Post-process a GAM containing non-standard effects
#' @name postproc_gam_nl
#' @rdname postproc_gam_nl
#' @export postproc_gam_nl
#'
postproc_gam_nl <- function(o, info) {
  
  ne <- length(info$type)
  
  for (ii in seq_len(ne)) {
    types <- info$type[[ii]]
    
    # Coefficients of nested smooth must be saved in smooth object
    if ("nested" %in% types) {
      
      ism <- info$extra[[ii]]$ism
      sii <- o$smooth[[ism]]
      si  <- sii$xt$si
      alpha_len <- length(si$alpha)
      si$alpha <- coef(o)[info$iec[[ii]][seq_len(alpha_len)]]
      
      # special case si_nexpsm, update alpha with positive constrain
      has_si_nexpsm <- "si_nexpsm" %in% types
      if (has_si_nexpsm) {
        n_nexp <- si$n_nexp
        n_si   <- si$n_si
        si$alpha_nexp <- as.vector(si$alpha[seq_len(n_nexp)])
        si$alpha_si   <- as.vector(si$alpha[(n_nexp + 1):(n_si + n_nexp)])
        
        if (isTRUE(si$positive_si)) {
          si$alpha_si_inner <- si$alpha_si             # optimized alpha_si
          si$alpha_si_true  <- exp(si$alpha_si)        # positive alpha_si
        }
      }
      sii$xt$si <- si 
      
      # Inner smooth must be centered using original data
      needs_base_xm <- is.null(si$xm)
      if (needs_base_xm || has_si_nexpsm) {
        inner_mean <- mean(attr(Predict.matrix.nested(sii, data = o$model), "inner_linpred_unscaled"))
        if (needs_base_xm) {
          si$xm <- inner_mean
        }
        if (has_si_nexpsm) {
          si$xm <- c(si$xm, nexp = inner_mean)
        }
        sii$xt$si$xm <- si$xm 
      }
      
      jacobian <- get_jacobian.nested(sii, data = o$model, param = coef(o)[info$iec[[ii]]])
      sii$xt$jacobian <- jacobian$JJ
      sii$xt$z        <- jacobian$z   # data after linear trans before exp smooth, useful for plot
      sii$xt$xsm_raw  <- jacobian$xsm_raw  # smoothed s_t, useful for plot(without center)
      sii$xt$xa       <- jacobian$xa # smoothed s_t, useful for plot
      o$smooth[[ism]] <- sii
    }
  }
  
  o$linear.predictors <- predict.gamnl(o, type = "link")
  o$fitted.values     <- predict.gamnl(o, type = "response")
  
  return(o)
}