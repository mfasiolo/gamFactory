#'
#' Function for setting up link functions
#'
#' @description Extension of \code{stats::make.link} (used internally by \code{mgcv})
#'              providing several additional link functions, needed to constrain GAM
#'              parameters to a bounded interval or half-line. If \code{link} is a name
#'              recognised by \code{stats::make.link} (e.g. \code{"identity"},
#'              \code{"log"}, \code{"logit"}, ...) that link is returned as usual;
#'              otherwise \code{link} is parsed as one of the custom links below, each
#'              specified as a string \code{"name(params)"}, e.g. \code{"logitab(-1, 2)"}.
#'
#' @param link a character string indicating the link function to be used: either a name
#'             recognised by \code{stats::make.link}, or one of the custom links below.
#'
#' @details The custom links currently available, with \eqn{\mu} the parameter on its
#'          natural scale and \eqn{\eta} the unconstrained linear predictor scale, are:
#' \describe{
#'   \item{\code{logitab(a, b)}}{Constrains \eqn{\mu \in (a, b)}, via a re-scaled logit:
#'         \deqn{\eta = \mathrm{logit}\left(\frac{\mu - a}{b - a}\right), \qquad
#'               \mu = a + (b - a) \, \mathrm{logit}^{-1}(\eta).}
#'         Used, e.g., for the shape parameter of the GEV distribution in
#'         \link{bundle_gev} (\code{"logitab(-1, 0.5)"}).}
#'   \item{\code{loginva(a)}}{Constrains \eqn{1/\mu \in (a, \infty)}, i.e.
#'         \eqn{\mu \in (0, 1/a)}:
#'         \deqn{\eta = \log(1/\mu - a), \qquad \mu = \frac{1}{\exp(\eta) + a}.}
#'         Used, e.g., for the inverse-scale parameter in \link{bundle_gaussian}
#'         (\code{"loginva(0.01)"}).}
#'   \item{\code{logea(a)}}{Constrains \eqn{\exp(\mu) \in (a, \infty)}:
#'         \deqn{\eta = \log(\exp(\mu) - a), \qquad \mu = \log(\exp(\eta) + a).}
#'         Used, e.g., for the kurtosis parameter in \link{bundle_shash}
#'         (\code{"logea(0.01)"}).}
#'   \item{\code{logea2(loga)}}{The same link as \code{logea(a)}, but parametrised by
#'         \code{loga = log(a)} instead of \code{a} itself, and evaluated in a
#'         numerically stable way. Constrains \eqn{\mu \in (loga, \infty)}, via a shifted
#'         softplus:
#'         \deqn{\mu = loga + \log(1 + \exp(\eta - loga)), \qquad
#'               \eta = loga + \log(\exp(\mu - loga) - 1).}
#'         Since \eqn{loga + \log(1 + \exp(\eta - loga)) = \log(\exp(loga) + \exp(\eta))
#'         = \log(a + \exp(\eta))}, this is mathematically the same link as
#'         \code{logea(exp(loga))}; the two agree to machine precision wherever both are
#'         finite (see examples). \code{logea2} exists because it evaluates
#'         \eqn{\exp(\eta - loga)} instead of \eqn{\exp(\eta)} directly, so unlike
#'         \code{logea} it does not overflow to \code{Inf} for large \eqn{\eta} or
#'         \eqn{\mu}. Used by default for the log-dispersion parameters in
#'         \link{bundle_gammals} and \link{bundle_gumbls} (\code{"logea2(-7)"}, i.e.
#'         \code{a = exp(-7)}).}
#' }
#'
#' @return An object of class \code{"link-glm"}, see \link[stats]{make.link}.
#' @name make_link
#' @rdname make_link
#' @importFrom gsubfn strapplyc
#' @export make_link
#' @examples
#' library(gamFactory)
#'
#' # logitab: mu constrained to (-1, 2)
#' lk <- make_link("logitab(-1, 2)")
#' range(lk$linkinv(seq(-50, 50, length.out = 200))) # always within (-1, 2)
#'
#' # loginva: 1/mu constrained to (2, Inf), i.e. mu in (0, 0.5)
#' lk <- make_link("loginva(2)")
#' range(lk$linkinv(seq(-50, 50, length.out = 200))) # always within (0, 0.5)
#'
#' # logea: exp(mu) constrained to (0.01, Inf)
#' lk <- make_link("logea(0.01)")
#' mu <- lk$linkinv(seq(-50, 50, length.out = 200))
#' all(exp(mu) > 0.01)
#'
#' # logea2: numerically-stable version of logea(exp(loga))
#' lk_ea2 <- make_link("logea2(-7)")
#' lk_ea  <- make_link(paste0("logea(", exp(-7), ")"))
#' eta <- seq(-5, 5, length.out = 100)
#' max(abs(lk_ea2$linkinv(eta) - lk_ea$linkinv(eta))) # ~ 0: same link
#' lk_ea2$linkinv(1000) # finite
#' lk_ea$linkinv(1000)  # overflows to Inf
#'
make_link <- function(link){
  
  # If link is in mgcv this will work ...
  out <- tryCatch(make.link(link), error = function(e) e)
  
  # Otherwise it is defined directly here
  if( inherits(out, "error") ){
    
    # NOTE: "logea2" contains "logea" as a substring, so it must be matched before
    # "logea" below (hence the if/else if chain, ordered most-specific-first).

    # Get limits (a, b) of logitab link
    if( grepl("logitab", link, fixed=TRUE) ){

      tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
      l1 <- tmp[1]
      l2 <- tmp[2]
      link <- "logitab"

    } else if( grepl("loginva", link, fixed=TRUE) ){
      # Get limit (a, \infinity) on 1/mu of loginva link

      tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
      l1 <- tmp[1]
      link <- "loginva"

    } else if( grepl("logea2", link, fixed=TRUE) ){
      # Get log-lower-limit loga of logea2 link (loga = log(a): same link as logea(a),
      # but parametrised by loga instead of a; used by default for the scale parameters
      # of mgcv::gammals and mgcv::gumbls). We extract the number from inside the
      # parentheses only, since the numeric-extraction regex below would otherwise also
      # match the "2" in "logea2" itself.
      inside <- sub(".*\\((.*)\\).*", "\\1", link)
      tmp <- as.numeric(strapplyc(inside, "[-+.e0-9]*\\d")[[1]])
      l1 <- tmp[1]
      link <- "logea2"

    } else if( grepl("logea", link, fixed=TRUE) ){
      # Get limit (a, \infinity) on exp(mu) of logea link

      tmp <- sort( as.numeric(strapplyc(link, "[-+.e0-9]*\\d")[[1]]) )
      l1 <- tmp[1]
      link <- "logea"

    }

    switch(link,
           logitab = {
             linkfun <- eval(parse(text=paste("function(mu) binomial()$linkfun((mu - ", 
                                              l1, ") / ", l2 - l1, ")", sep='')))
             mu.eta <- eval(parse(text=paste("function(eta) binomial()$mu.eta(eta) * ", l2 - l1, sep=''))) 
             linkinv <- eval(parse(text=paste("function(eta) binomial()$linkinv(eta) * ",
                                              l2 - l1, " + ", l1, sep='')))
             valideta <- function(eta) all( is.finite(eta) )
           }, 
           loginva = {
             linkfun <-  eval(parse(text=paste("function(mu) log(1/mu -",l1,")")))
             mu.eta <- eval(parse(text=paste("function(eta) { ee <- exp(eta); -ee/(ee +",l1,")^2 }")))
             linkinv <- eval(parse(text=paste("function(eta) 1/(exp(eta) +",l1,")")))
             valideta <- function(eta) all( is.finite(eta) )
           }, 
           logea = {
             linkfun <-  eval(parse(text=paste("function(mu) log(exp(mu) - ",l1,")", sep='')))
             mu.eta <-  eval(parse(text=paste("function(eta) { ee <- exp(eta); ee/(ee +",l1,") }")))
             linkinv <- eval(parse(text=paste("function(eta) log(exp(eta) +",l1,")", sep='')))
             valideta <- function(eta) all( is.finite(eta) )
           },
           logea2 = {
             # Numerically-stable reformulation of logea(a) with a = exp(l1), i.e. l1 is
             # loga = log(a). Mathematically linkinv(eta) = log(exp(eta) + exp(l1)), but
             # evaluated as l1 + log1p(exp(eta - l1)) so that it never overflows for
             # large eta or mu (unlike logea, which computes exp(mu)/exp(eta) directly).
             linkfun <- eval(parse(text=paste0(
               "function(mu) {\n",
               "  eta <- mub <- mu - (", l1, ")\n",
               "  ii <- mub < .Machine$double.eps\n",
               "  if (any(ii)) eta[ii] <- log(.Machine$double.eps) + (", l1, ")\n",
               "  jj <- mub > -log(.Machine$double.eps)\n",
               "  if (any(jj)) eta[jj] <- mub[jj] + (", l1, ")\n",
               "  jj <- !jj & !ii\n",
               "  if (any(jj)) eta[jj] <- log(exp(mub[jj]) - 1) + (", l1, ")\n",
               "  eta\n}")))
             mu.eta <- eval(parse(text=paste0(
               "function(eta) {\n",
               "  eta <- eta - (", l1, ")\n",
               "  ii <- eta < 0\n",
               "  eta <- exp(-eta * sign(eta))\n",
               "  if (any(ii)) { ei <- eta[ii]; eta[ii] <- ei / (1 + ei) }\n",
               "  ii <- !ii\n",
               "  if (any(ii)) eta[ii] <- 1 / (1 + eta[ii])\n",
               "  eta\n}")))
             linkinv <- eval(parse(text=paste0(
               "function(eta) {\n",
               "  mu <- eta\n",
               "  ii <- eta - (", l1, ") < -log(.Machine$double.eps)\n",
               "  if (any(ii)) mu[ii] <- (", l1, ") + log1p(exp(eta[ii] - (", l1, ")))\n",
               "  mu\n}")))
             valideta <- function(eta) TRUE
           },
           stop(gettextf("%s link not recognised", sQuote(link)), domain = NA))
    
    environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- 
      environment(valideta) <- asNamespace("stats")
    
    out <- structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
                          valideta = valideta, name = link), class = "link-glm")
    
  }
  
  return( out )
}