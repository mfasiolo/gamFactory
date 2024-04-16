# We are trying to reuse the methods provided by 
# stats::poisson, stats::binomial, ...
.resid_exp_fam <- function(object, type, fam){
  y <- drop(object$y)
  mu <- object$fitted.values
  Ey <- mu
  # Expected value of y is mu, unless another expected value function is provided
  if(!is.null(fam$exp_val)){
   Ey <- fam$exp_val(mu)
  }
  wts <- object$prior.weights
  if(type == "deviance"){
    rsd <- fam$dev.resid(y, mu, wts)
    s <- attr(rsd, "sign")
    if (is.null(s)){
      s <- sign(y - Ey)
    }
    rsd <- sqrt(pmax(rsd, 0)) * s
  } else {
    rsd <- y - Ey
    if(type == "pearson"){
      rsd <- rsd * sqrt(wts / fam$var(mu))
    }
  }
  return(rsd)
}