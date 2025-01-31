#' Processing the summary to be printed
#' @description  This function allows to modify the summary.gam function, to adjust p-values for nested effects
#'
#' @param object an object of class gamnl
#' @param intercept set to TRUE if you want print the intercept
#' @param print print the summary
#' @param ... further arguments to be passed to gam
#'
#' @return The summary is printed
#' @export summary.gamnl
#' @export
#'
#' @importFrom mgcv summary.gam
#' @importFrom mgcv print.summary.gam
#' @name summary.gamnl
#' @rdname summary
summary.gamnl <- function(object, intercept = FALSE, print = TRUE, ...){
  
  sgam <- summary.gam(object)
  
  sgam.st <- sgam$s.table
  lab_st <- rownames(sgam.st)
  terms_nl <- unlist(sapply(1:length(object$smooth), 
                    function(ii) { 
                      ss <- object$smooth[[ii]]
                      term <- ifelse("nested" %in% class(ss), ii, 0)
                      names(term) <- ss$label
                      return(term)
                      }))
  for (ii in 1:length(lab_st)) {
    if(terms_nl[lab_st[ii]]>0){
      sm <- object$smooth[[terms_nl[lab_st[ii]]]]
      ss <- predict(object, type = "terms", terms = lab_st[ii])
      Jac <- sm$xt$jacobian
      if(sm$rank[1] == sm$bs.dim - length(sm$xt$si$alpha)){ # if full rank use Ve
        V <- object$Ve[sm$first.para:sm$last.para,sm$first.para:sm$last.para]
      } else{
        V <- object$Vp[sm$first.para:sm$last.para,sm$first.para:sm$last.para]
      }
      pval <- .sumChi2(ss, V, Jac)
      
      sgam.st[ii, "Chi.sq"] <- pval$stat
      sgam.st[ii, "p-value"] <- pval$pval
      lab_st[ii] <- gsub("s\\(", "s_nest\\(", lab_st[ii])
    }
  }
  sgam$chi.sq <- sgam.st[,"Chi.sq"]
  sgam$s.pv <- sgam.st[,"Chi.sq"]
  rownames(sgam.st) <- names(sgam$chi.sq) <- lab_st
  sgam$s.table <- sgam.st
  
  if(print){
    print.summary.gam(sgam)
  }
  ret <- with(sgam, list(p.coeff = p.coeff, se = se, p.t = p.t, p.pv = p.pv,
                         residual.df = residual.df, m = m, chi.sq = chi.sq, s.pv = s.pv,
                         scale = dispersion, r.sq = r.sq, family = object$family,
                         formula = object$formula, n = nobs, dev.expl = dev.expl,
                         edf = edf, dispersion = dispersion, pTerms.pv = pTerms.pv,
                         pTerms.chi.sq = pTerms.chi.sq, pTerms.df = pTerms.df,
                         cov.unscaled = cov.unscaled, cov.scaled = cov.scaled,
                         p.table = p.table, pTerms.table = pTerms.table, s.table = sgam.st,
                         method = object$method, sp.criterion = object$gcv.ubre,
                         rank = object$rank, np = length(object$coefficients)))
  class(ret) <- "summary.gamnl"
  return(invisible(ret))
}


