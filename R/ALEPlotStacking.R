#' Accumulated Local Effects (ALE) Plots of a 
#' stacking of predictive distribution model.
#'
#' @param X Data frame with the predictor variables used in the formula, 
#' not including the response.
#' @param X.model The fitted model.
#' @param J A numeric scalar indicating the column of the predictor 
#' in the model matrix for which the ALE plot will be calculated. 
#' At the moment only effects of single covariates are implemented (no interactions)
#' @param K Number of intervals into which the predictor range is divided 
#' when calculating the ALE plot effects
#' @param which.alpha A numeric scalar indicating the category for which the 
#' effect will be plotted.
#' @param confidence.intervals A logical value. If TRUE, 
#' approximate posterior credible intervals based on the Delta method are plotted.
#' Intervals are available only for numeric variables, not for factor variables.
#' @return The ALE plot of the covariate effects over each category/expert.
#' @export ALEPlotStacking
#' @name ALEPlotStacking
#' @rdname ALEPlotStacking
#' 
ALEPlotStacking <- function(X, X.model, J, which.alpha, K = 40,
                            plot_output = TRUE,
                            confidence.intervals = TRUE,
                            ylim = NULL) {
  
  pred.fun <- function(X.model, newdata) {
    nuHat <- predict(X.model, newdata = newdata)
    aHat <- cbind(1, exp(nuHat)) / rowSums(cbind(1, exp(nuHat)))
    aHat
  }
  
  N = dim(X)[1]
  d = dim(X)[2]
  
  if (class(X[, J]) == "factor") {
    X[, J] <- droplevels(X[, J])
    x.count <- as.numeric(table(X[, J]))
    x.prob <- x.count/sum(x.count)
    K <- nlevels(X[, J])
    D.cum <- matrix(0, K, K)
    D <- matrix(0, K, K)
    for (j in setdiff(1:d, J)) {
      if (class(X[, j]) == "factor") {
        A = table(X[, J], X[, j])
        A = A/x.count
        for (i in 1:(K - 1)) {
          for (k in (i + 1):K) {
            D[i, k] = sum(abs(A[i, ] - A[k, ]))/2
            D[k, i] = D[i, k]
          }
        }
        D.cum <- D.cum + D
      }
      else {
        q.x.all <- quantile(X[, j], probs = seq(0, 
                                                1, length.out = 100), 
                            na.rm = TRUE, names = FALSE)
        x.ecdf = tapply(X[, j], X[, J], ecdf)
        for (i in 1:(K - 1)) {
          for (k in (i + 1):K) {
            D[i, k] = max(abs(x.ecdf[[i]](q.x.all) - 
                                x.ecdf[[k]](q.x.all)))
            D[k, i] = D[i, k]
          }
        }
        D.cum <- D.cum + D
      }
    }
    D1D <- cmdscale(D.cum, k = 1)
    ind.ord <- sort(D1D, index.return = T)$ix
    ord.ind <- sort(ind.ord, index.return = T)$ix
    levs.orig <- levels(X[, J])
    levs.ord <- levs.orig[ind.ord]
    x.ord <- ord.ind[as.numeric(X[, J])]
    row.ind.plus <- (1:N)[x.ord < K]
    row.ind.neg <- (1:N)[x.ord > 1]
    X.plus <- X
    X.neg <- X
    X.plus[row.ind.plus, J] <- levs.ord[x.ord[row.ind.plus] + 
                                          1]
    X.neg[row.ind.neg, J] <- levs.ord[x.ord[row.ind.neg] - 
                                        1]
    y.hat <- pred.fun(X.model = X.model, newdata = X)
    y.hat.plus <- pred.fun(X.model = X.model, newdata = X.plus[row.ind.plus, 
                                                               ])
    y.hat.neg <- pred.fun(X.model = X.model, newdata = X.neg[row.ind.neg, 
                                                             ])
    Delta.plus <- y.hat.plus[, which.alpha] - y.hat[row.ind.plus, which.alpha]
    Delta.neg <- y.hat[row.ind.neg, which.alpha] - y.hat.neg[, which.alpha]
    Delta <- as.numeric(tapply(c(Delta.plus, Delta.neg), 
                               c(x.ord[row.ind.plus], x.ord[row.ind.neg] - 
                                   1), mean))
    fJ <- c(0, cumsum(Delta))
    fJ = fJ - sum(fJ * x.prob[ind.ord])
    x <- levs.ord
    
    p <- ggplot(data.frame(x, fJ)) +
      theme_bw() +
      geom_col(aes(x = x, y = fJ)) +
      xlab(paste("x_", J, " (", names(X)[J], ")", sep = "")) +
      ylab(paste("f_", J, "(x_", J, ")", sep = "")) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5))
    
    if (!is.null(ylim)) {
      p <- p + ylim(ylim)
    }
    
    # barplot(fJ, names = x, xlab = paste("x_", J, " (", 
    #                                     names(X)[J], ")", sep = ""), 
    #         ylab = paste("f_", J, "(x_", J, ")", sep = ""), las = 3)
    se <- NULL ## confidence intervals not yet implemented
  }
  else if (class(X[, J]) == "numeric" | class(X[, J]) == "integer") {
    z = c(min(X[, J]), 
          as.numeric(quantile(X[, J], seq(1/K, 1, length.out = K), type = 1)))
    z = unique(z)
    K = length(z) - 1
    fJ = numeric(K)
    a1 = as.numeric(cut(X[, J], breaks = z, include.lowest = TRUE))
    X1 = X
    X2 = X
    X1[, J] = z[a1]
    X2[, J] = z[a1 + 1]
    y.hat1 = pred.fun(X.model = X.model, newdata = X1)
    y.hat2 = pred.fun(X.model = X.model, newdata = X2)
    
    an1 <- sapply(1:(ncol(y.hat1) - 1), function(hh) {
      y.hat1[, which.alpha] * (as.numeric(which.alpha == hh + 1) -
                                 y.hat1[, hh + 1])
    })
    an2 <- sapply(1:(ncol(y.hat2) - 1), function(hh) {
      y.hat2[, which.alpha] * (as.numeric(which.alpha == hh + 1) -
                                 y.hat2[, hh + 1])
    })
    
    Z1 <- model.matrix(X.model, newdata = X1)
    Z2 <- model.matrix(X.model, newdata = X2)
    lpi <- attr(Z1, "lpi")
    
    J1 <- lapply(seq_along(lpi), function(ii) {
      Z1[, lpi[[ii]]] * an1[, ii]
    }) %>% do.call(cbind, .)
    
    J2 <- lapply(seq_along(lpi), function(ii) {
      Z2[, lpi[[ii]]] * an2[, ii]
    }) %>% do.call(cbind, .)
    
    # J1 <- Z1 * an1
    # J2 <- Z2 * an2
    
    Fb <- J2 - J1
    Fb <- apply(Fb, 2, function(fb) as.numeric(tapply(fb, a1, mean)))
    Fb <- apply(Fb, 2, cumsum)
    b1 <- as.numeric(table(a1))
    
    ## ALEPlot centers the effects differently and adds 0 to the cumsum
    
    D <- apply(Fb, 2, function(fb) fb - sum(fb * b1) / sum(b1))
    
    Delta = y.hat2[, which.alpha] - y.hat1[, which.alpha]
    Delta = as.numeric(tapply(Delta, a1, mean))
    
    # fJ = c(0, cumsum(Delta))
    fJ = cumsum(Delta)
    
    # fJ = fJ - sum((fJ[1:K] + fJ[2:(K + 1)])/2 * b1)/sum(b1)
    fJ = fJ - sum(fJ * b1) / sum(b1)
    x <- z
    
    x <- x[- 1]
    
    covF <- D %*% X.model$Vp %*% t(D)
    se <- sqrt(diag(covF))
    
    if (plot_output) {
      p <- ggplot(data.frame(x, fJ)) +
        theme_bw() +
        geom_line(aes(x = x, y = fJ)) +
        xlab(paste("x_", J, " (", names(X)[J], ")", sep = "")) +
        ylab(paste("f_", J, "(x_", J, ")", sep = ""))
      if (!is.null(ylim)) {
        p <- p + ylim(ylim)
      } else {
        p <- p +
          ylim(c(min(fJ - 1.96 * se), max(fJ + 1.96 * se)))
      }
      
      
      # plot(x, fJ, type = "l",
      #      xlab = paste("x_", J, " (", names(X)[J], ")", sep = ""),
      #      ylab = paste("f_", J, "(x_", J, ")", sep = ""),
      #      
      #      ylim = c(min(fJ - 1.96 * se), max(fJ + 1.96 * se))
      # )
      
      if (confidence.intervals) {
        p <- p + 
          geom_line(aes(x, fJ + 1.96 * se), lty = 2) +
          geom_line(aes(x, fJ - 1.96 * se), lty = 2)
        
        # lines(x, fJ + 1.96 * se, lty = 2)
        # lines(x, fJ - 1.96 * se, lty = 2)
      }
      p
    }
  }
  else print("error:  
             class(X[,J]) must be either factor or numeric or integer")
  
  list(K = K, x.values = x, f.values = fJ, se = se, p = p)
}