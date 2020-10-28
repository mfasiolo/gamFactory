#'
#' Build linear predictors
#' 
#' @param eff list of effects, see e.g. [eff_stand] or [eff_si].
#' @param iel vectors indicating to which linear predictor each effect belongs to.
#' @param iec list of vectors indicating which regression coefficients belong to each effect.
#' 
#' @rdname linpreds
#' @export linpreds
#'
linpreds <- function(eff, iel, iec){
  
  # Number of linear predictors
  nlp <- length( table(iel) )

  # List of vectors indicating which components belong to each linear predictor
  ile <- lapply(1:nlp, function(ii) which(iel == ii))
  
  # Number of components
  nc <- length( eff )
  
  # Evaluate all the linear predictors
  eval <- function(param, deriv = 0){
    
    # Looping over LPs, then components of that LP
    eff <- lapply(1:nc, 
                  function(ii){
                    eff[[ii]]$eval(param = param[iec[[ii]]], deriv = deriv)
                  })
    
    o <- linpreds(eff = eff, iel = iel, iec = iec)
    
    # For each linear predictor, sum its evaluated effects 
    o$f <- lapply(ile, 
                    function(ii){
                      Reduce("+", lapply(ii, function(kk) eff[[kk]]$f ))
                    })
    
    o$param <- param
    o$deriv <- deriv
    
    return( o )
    
  }
  
  out <- structure(list("eff" = eff, "eval" = eval, "nc" = nc, "nlp" = nlp, "iel" = iel, "iec" = iec, "ile" = ile), 
                   class = "multiLP")
  
  return( out )
  
}