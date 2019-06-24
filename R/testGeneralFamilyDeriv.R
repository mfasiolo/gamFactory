#' Test derivatives of a general family to be used in mgcv
#'
#' @param fam A general family
#' @param nsim Number of parameters to simulate 
#' (both for regression coefficients and smoothing parameters) 
#' at which calculate derivatives
#' @param ndata Number of data to generate
#' @param p Number of regression coefficients for each linear predictor
#' @param seed Set a seed for reproducibility
#' @param ntheta Number of additional parameters of the family (to be implemented)
#'
#' @return A comparison between derivatives calculated according to the ll function of the family
#' and numerical derivatives calculated using jacobian function
#' @export
#'
testGeneralFamilyDeriv <- function(fam,
                            nsim = 30,
                            ndata = 1e2,
                            p = 2,
                            seed = 0,
                            nlp = 3,
 ## add code to manage additional parameters theta depending on the family.
 ## The problem is that theta parameters are added to the regression coefficients
 ## then this function is not able to manage the correct number of 
 ## regression coefficients, linear predictors, response variables etc.
 
 ## this tests the same formula for each linear predictor both for
 ## stacking of means and stacking of predictive distributions.
 ## different formulas per predictor not yet implemented
                            ntheta = NULL 
                            ) {
  
  ## simulate data
  set.seed(seed)
  y <- matrix(rnorm(ndata * nlp), nrow = ndata, ncol = nlp)
  x <- matrix(rnorm(ndata * p), nrow = ndata, ncol = p * nlp)
  attr(x, "lpi") <- lapply(1:nlp, function(ii) ((ii - 1) * p + 1):(ii * p))
  npars <- p * nlp
  
  ## specific for stackPredictiveFamily
  if (identical(fam, stackPredictiveFamily)) {
    K <- nlp + 1
    X <- matrix(rnorm(ndata * K), nrow = ndata, ncol = K)
    fam <- stackPredictiveFamily(X)
  }
    
  ## specific for stackFamily
  if (identical(fam, stackFamily)) {
    K <- nlp
    X <- matrix(rnorm(ndata * K), nrow = ndata, ncol = K)
    fam <- stackFamily(X, familyDeriv = createGaussian(y[, 1]))
    npars <- npars + 1 ## gaussian family has one additional parameter
  }
    
  if (class(fam)[1] == "function") fam <- fam()
  if (!("general.family" %in% class(fam))) 
    stop("Only general families can be tested with this function")
  
  ## first and second derivatives wrt regression coefficients
  dl <- function(coef) {
    fam$ll(y, x, coef, deriv = 1, family = fam)$l
  }
  
  dlb <- function(coef) {
    fam$ll(y, x, coef, deriv = 1, family = fam)$lb
  }
  
  dlbNum <- function(coef) {
    jacobian(dl, coef)
  }
  
  dlbb <- function(coef) {
    fam$ll(y, x, coef, deriv = 1, family = fam)$lbb %>%
      as.numeric
  }
  
  dlbbNum <- function(coef) {
    jacobian(dlb, coef) %>% as.numeric
  }
  
  beta <- matrix(rnorm(nsim * npars), nrow = nsim, ncol = npars)
  
  ## Solves the problem for gaulss family,
  ## where the number of linear predictors does not equal
  ## the number of columns of y
  if (fam$family == "gaulss") y <- y[, 1]
  
  funsBeta <- c("dlb", "dlbNum", "dlbb", "dlbbNum")
  derivListBeta <-
    lapply(funsBeta, function(fun) {
      apply(beta, 1, function (b) do.call(fun, 
                                          list(coef = b))
            )
    })
  names(derivListBeta) <- funsBeta
  

  ## Derivatives of Hessian wrt smoothing parameters
  # derivative of the hessian of the log-likelihood wrt rho
  
  pen_ll <- function(beta, lam) {
    
    penalty <- lapply(1:nsp, function(ii) {
      lam[ii] * t(beta) %*% S[[ii]] %*% beta
    }) %>%
      Reduce("+", .)
    
    pen_gr <- lapply(1:nsp, function(ii) {
      lam[ii] * S[[ii]] %*% beta
    }) %>%
      Reduce("+", .)
    
    pen_hes <- lapply(1:nsp, function(ii) {
      lam[ii] * S[[ii]]
    }) %>%
      Reduce("+", .)
    
    derObj <- fam$ll(y, x, beta, family = fam, deriv = 1)
    
    out <- - derObj$l + .5 * penalty
    attr(out, "gradient") <- - derObj$lb + pen_gr
    attr(out, "hessian") <- - derObj$lbb + pen_hes
    out
  }
  
  beta0 <- rep(0, npars)
  nsp <- nlp
  nnn <- as.integer(npars / nsp)
  S <- lapply(1:nsp, function (ii) {
    S_ii <- matrix(0, nrow = npars, ncol = npars)
    diag(S_ii)[1:nnn + nnn * (ii - 1)] <- 1
    S_ii
  })
  
  # function optimizing beta given the smoothing parameter
  d3r <- function(r) {

    lam <- exp(r)
    # Inner step: optimize pen l wrt beta
    
    optimization <- nlm(pen_ll, beta0, lam = lam, hessian = T 
                        # gradtol = 1e-10
                        )
    
    betaHat <- optimization$estimate # beta hat as function of rho
    negHessian <- optimization$hessian
    invNegH <- solve(optimization$hessian)
    d1br <- lapply(1:nsp, function(ii){
      - lam[ii] * invNegH %*% S[[ii]] %*% betaHat
    }) %>% do.call("cbind", .)
    list(betaHat = betaHat, d1br = d1br, negHessian = negHessian)
  }
  
  d1h <- function(sp) {
    
    d3rObj <- d3r(sp)
    betaHat <- d3rObj$betaHat
    d1b <- d3rObj$d1br

    fam$ll(y, x, betaHat, family = fam, d1b = d1b, deriv = 3)$d1H %>% 
      setNames(paste0("sp", 1:nsp))
  }
  
  dd1HNum <- function(sp) {
    jacobian(function(x) dlbb(d3r(x)$betaHat), sp) %>% 
      as.data.frame %>% 
      lapply(function(x) matrix(x, nrow = nlp * p, ncol = nlp * p)) %>% 
      setNames(paste0("sp", 1:nsp))
  }
  
  nsp <- nlp
  
  organizeD1H <- function(d1H) {
    lapply(1:nsim, function(ii) {
      lapply(1:nsp, function(jj) {
        data.frame(d1H = as.numeric(d1H[[ii]][[jj]]),
                   isDiag = as.logical(diag(p * nsp)),
                   nsim = ii,
                   nsp = jj
        )
      }) %>% 
        bind_rows
    }) %>% 
      bind_rows %>% 
      mutate(nsim = factor(nsim),
             nsp = factor(nsp))
  }
  
  r <- matrix(rnorm(nsim * nsp), nrow = nsim, ncol = nsp)
  d1HEx <- apply(r, 1, d1h) %>% 
    setNames(paste0("sim", 1:nsim)) %>% 
    organizeD1H
  
  d1HNum <- apply(r, 1, dd1HNum) %>% 
    setNames(paste0("sim", 1:nsim)) %>% 
    organizeD1H
  
  derivSp <- data.frame(d1HEx, Numeric = d1HNum$d1H) %>% 
    rename(Exact = d1H) %>% 
    mutate(d1HDiff = Exact - Numeric,
           relDiff = abs(Exact - Numeric) / 
             apply(abs(data.frame(Exact, Numeric)), 1, max))
  
  columns <- list(c("dlb", "dlbNum"), c("dlbb", "dlbbNum"))
  names(columns) <- 
    c("First derivatives wrt beta", 
      "Second derivatives wrt beta")
  lapply(seq_along(columns), function(i, col) {
    data.frame(Exact = derivListBeta[[columns[[i]][1]]] %>% as.numeric, 
               Numeric = derivListBeta[[columns[[i]][2]]] %>% as.numeric) %>%
      ggplot() +
      geom_point(aes(Exact, Numeric), size = 0.1) +
      geom_abline(col = "red") +
      ggtitle(col[[i]])
  }, 
  col = names(columns)
  ) %>% 
    do.call(grid.arrange, .)
  
  lapply(c("isDiag", "nsim", "nsp"), function(colu) {
    ggplot(derivSp) + 
      geom_point(aes_string("Exact", "Numeric", col = colu), size = 0.5) +
      geom_abline(col = "red") +
      ggtitle("First derivatives of Hessian wrt smoothing parameters")
  }) %>% 
    do.call(grid.arrange, .)

  
  relDiff <- function(x) abs(diff(x)) / max(abs(x))
  
  lb <- data.frame(EX = derivListBeta$dlb %>% as.numeric,
                    NUM = derivListBeta$dlbNum %>% as.numeric) %>% 
    mutate(DIFF = apply(., 1, relDiff))
  
  lbb <- data.frame(EX = derivListBeta$dlbb %>% as.numeric,
                    NUM = derivListBeta$dlbbNum %>% as.numeric) %>% 
    mutate(DIFF = apply(., 1, relDiff))
  
  d1H <- data.frame(EX = derivSp$Exact %>% as.numeric,
             NUM = derivSp$Numeric %>% as.numeric) %>%
    mutate(DIFF = apply(., 1, relDiff))

  list(lb = lb, lbb = lbb, d1H = d1H)
  
}

