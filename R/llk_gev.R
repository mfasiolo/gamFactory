##################
#'
#' Log-likelihood of the Generalized Extreme Value (GEV) distribution
#'
#' @description Log-likelihood of the location-scale-shape GEV distribution, and its
#'              derivatives with respect to \code{param = cbind(mu, rho, xi)}, where
#'              \code{mu} is the location, \code{rho = log(sigma)} is the log-scale and
#'              \code{xi} is the shape parameter. The algebra here is ported directly
#'              from the derivative system used internally by \code{mgcv::gevlss}, but
#'              re-expressed as derivatives w.r.t. \code{param} and returned in the list
#'              format used by [`gamFactory::llk_gaussian`] and friends.
#' @param y a vector of observations.
#' @param param a matrix (or list) with 3 columns (elements), containing the location,
#'              log-scale and shape parameters, in this order.
#' @param deriv integer between 0 and 4 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-4
#'              additionally return \code{d1}-\code{d4}.
#' @rdname llk_gev
#' @export llk_gev
#' @examples
#' library(gamFactory)
#'
#' n <- 10
#' param <- cbind(mu = rnorm(n), rho = rnorm(n), xi = runif(n, -0.3, 0.3))
#' y <- param[ , 1] + exp(param[ , 2]) * (((-log(runif(n)))^(-param[ , 3]) - 1) / param[ , 3])
#'
#' llk_gev(y = y, param = param, deriv = 4)
#'
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_gev(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_gev(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 1)$d1))
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_gev(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 2)$d2))
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_gev(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 3)$d3))
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_gev(y = y, param = matrix(param, n, 3, byrow = TRUE), deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = c(0, 0, 0.1), ord = 1:4)
#' @export
#'
llk_gev <- function(y, param, deriv = 0, ...) {

  if (is.list(param)) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 3) stop("Wrong number of parameters provided")

  n <- length(y)
  if (nrow(param) == 1 && n > 1) {
    param <- matrix(param, n, 3, byrow = TRUE)
  }

  mu  <- param[ , 1, drop = TRUE]
  rho <- param[ , 2, drop = TRUE]
  xi  <- param[ , 3, drop = TRUE]

        eps <- 1e-07
        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        exp1 <- exp(1)
        ymu <- y - mu
        aa0 <- (xi * ymu)/exp1^rho
        log.aa1 <- log1p(aa0)
        aa1 <- aa0 + 1
        aa2 <- 1/xi
        l0 <- (-aa2 * (1 + xi) * log.aa1) - 1/aa1^aa2 - rho

  out <- list("d0" = l0)

  if (deriv > 0) {
            bb1 <- 1/exp1^rho
            bb2 <- bb1 * xi * ymu + 1
            l1_1 <- (bb1 * (xi + 1))/bb2 - bb1 * bb2^((-1/xi) - 
                1)
            cc2 <- ymu
            cc0 <- bb1 * xi * cc2
            log.cc3 <- log1p(cc0)
            cc3 <- cc0 + 1
            l1_2 <- (-bb1 * cc2 * cc3^((-1/xi) - 1)) + (bb1 * 
                (xi + 1) * cc2)/cc3 - 1
            dd3 <- xi + 1
            dd6 <- 1/cc3
            dd7 <- log.cc3
            dd8 <- 1/xi^2
            l1_3 <- (-(dd8 * dd7 - bb1 * aa2 * cc2 * dd6)/cc3^aa2) + 
                dd8 * dd3 * dd7 - aa2 * dd7 - bb1 * aa2 * dd3 * 
                cc2 * dd6

    out[["d1"]] <- list(l1_1, l1_2, l1_3)

    if (deriv > 1) {
            ee1 <- 1/exp1^(2 * rho)
            ee3 <- -1/xi
            l2_1 <- ee1 * (ee3 - 1) * xi * aa1^(ee3 - 2) + 
                (ee1 * xi * (xi + 1))/aa1^2
            ff7 <- ee3 - 1
            l2_2 <- bb1 * cc3^ff7 + ee1 * ff7 * xi * cc2 * 
                cc3^(ee3 - 2) - (bb1 * dd3)/cc3 + (ee1 * xi * 
                dd3 * cc2)/cc3^2
            gg7 <- -aa2
            l2_3 <- (-bb1 * cc3^(gg7 - 1) * (log.cc3/xi^2 - 
                bb1 * aa2 * cc2 * dd6)) + ee1 * cc2 * cc3^(gg7 - 
                2) + bb1 * dd6 - (ee1 * (xi + 1) * cc2)/cc3^2
            hh4 <- cc2^2
            l2_4 <- bb1 * cc2 * cc3^ff7 + ee1 * ff7 * xi * 
                hh4 * cc3^(ee3 - 2) - (bb1 * dd3 * cc2)/cc3 + 
                (ee1 * xi * dd3 * hh4)/cc3^2
            l2_5 <- (-bb1 * cc2 * cc3^(gg7 - 1) * (log.cc3/xi^2 - 
                bb1 * aa2 * cc2 * dd6)) + ee1 * hh4 * cc3^(gg7 - 
                2) + bb1 * cc2 * dd6 - (ee1 * (xi + 1) * hh4)/cc3^2
            jj08 <- 1/cc3^2
            jj12 <- 1/xi^3
            jj13 <- 1/cc3^aa2
            l2_6 <- (-jj13 * (dd8 * dd7 - bb1 * aa2 * cc2 * 
                dd6)^2) - jj13 * (ee1 * aa2 * hh4 * jj08 + 2 * 
                bb1 * dd8 * cc2 * dd6 - 2 * jj12 * dd7) - 2 * 
                jj12 * dd3 * dd7 + 2 * dd8 * dd7 + 2 * bb1 * 
                dd8 * dd3 * cc2 * dd6 - 2 * bb1 * aa2 * cc2 * 
                dd6 + ee1 * aa2 * dd3 * hh4 * jj08

      out[["d2"]] <- list(l2_1, l2_2, l2_3, l2_4, l2_5, l2_6)

      if (deriv > 2) {
            kk1 <- 1/exp1^(3 * rho)
            kk2 <- xi^2
            l3_1 <- (2 * kk1 * kk2 * (xi + 1))/aa1^3 - kk1 * 
                (ee3 - 2) * (ee3 - 1) * kk2 * aa1^(ee3 - 3)
            ll5 <- (xi * cc2)/exp1^rho + 1
            ll8 <- ee3 - 2
            l3_2 <- (-2 * ee1 * ff7 * xi * ll5^ll8) - kk1 * 
                ll8 * ff7 * kk2 * cc2 * ll5^(ee3 - 3) - (2 * 
                ee1 * xi * dd3)/ll5^2 + (2 * kk1 * kk2 * dd3 * 
                cc2)/ll5^3
            mm10 <- cc3^(gg7 - 3)
            mm11 <- gg7 - 2
            mm12 <- cc3^mm11
            l3_3 <- ee1 * (gg7 - 1) * xi * mm12 * (log.cc3/xi^2 - 
                (bb1 * aa2 * cc2)/cc3) - ee1 * mm12 - kk1 * mm11 * 
                xi * cc2 * mm10 + kk1 * cc2 * mm10 + ee1 * dd3 * 
                jj08 + ee1 * xi * jj08 - (2 * kk1 * xi * dd3 * 
                cc2)/cc3^3
            l3_4 <- (-bb1 * cc3^ff7) - 3 * ee1 * ff7 * xi * 
                cc2 * cc3^ll8 - kk1 * ll8 * ff7 * kk2 * hh4 * 
                cc3^(ee3 - 3) + (bb1 * dd3)/cc3 - (3 * ee1 * 
                xi * dd3 * cc2)/cc3^2 + (2 * kk1 * kk2 * dd3 * 
                hh4)/cc3^3
            oo10 <- gg7 - 1
            oo13 <- log.cc3/xi^2
            l3_5 <- bb1 * cc3^oo10 * (bb1 * oo10 * cc2 * dd6 + 
                oo13) + ee1 * oo10 * xi * cc2 * mm12 * (bb1 * 
                mm11 * cc2 * dd6 + oo13) + ee1 * aa2 * cc2 * 
                mm12 + ee1 * oo10 * cc2 * mm12 - bb1 * dd6 + 
                2 * ee1 * dd3 * cc2 * jj08 + ee1 * xi * cc2 * 
                jj08 - (2 * xi * dd3 * cc2^2)/(exp1^(3 * rho) * 
                cc3^3)
            pp07 <- (-1/xi) - 1
            pp08 <- cc3^pp07
            l3_6 <- (-bb1 * pp08 * (bb1 * pp07 * cc2 * dd6 + 
                dd8 * dd7)^2) - bb1 * pp08 * ((-ee1 * pp07 * 
                hh4 * jj08) + 2 * bb1 * dd8 * cc2 * dd6 - (2 * 
                dd7)/xi^3) - 2 * ee1 * cc2 * jj08 + (2 * (xi + 
                1) * hh4)/(exp1^(3 * rho) * cc3^3)
            qq05 <- cc2^3
            l3_7 <- (-bb1 * cc2 * cc3^ff7) - 3 * ee1 * ff7 * 
                xi * hh4 * cc3^ll8 - kk1 * ll8 * ff7 * kk2 * 
                qq05 * cc3^(ee3 - 3) + (bb1 * dd3 * cc2)/cc3 - 
                (3 * ee1 * xi * dd3 * hh4)/cc3^2 + (2 * kk1 * 
                kk2 * dd3 * qq05)/cc3^3
            rr17 <- log.cc3/xi^2 - bb1 * aa2 * cc2 * dd6
            l3_8 <- bb1 * cc2 * cc3^oo10 * rr17 + ee1 * oo10 * 
                xi * hh4 * mm12 * rr17 - 2 * ee1 * hh4 * mm12 - 
                kk1 * mm11 * xi * qq05 * mm10 + kk1 * qq05 * 
                mm10 - bb1 * cc2 * dd6 + 2 * ee1 * dd3 * hh4 * 
                jj08 + ee1 * xi * hh4 * jj08 - (2 * kk1 * xi * 
                dd3 * qq05)/cc3^3
            l3_9 <- (-bb1 * cc2 * pp08 * (bb1 * pp07 * cc2 * 
                dd6 + dd8 * dd7)^2) - bb1 * cc2 * pp08 * ((-ee1 * 
                pp07 * hh4 * jj08) + 2 * bb1 * dd8 * cc2 * dd6 - 
                (2 * dd7)/xi^3) - 2 * ee1 * hh4 * jj08 + (2 * 
                (xi + 1) * cc2^3)/(exp1^(3 * rho) * cc3^3)
            tt08 <- 1/cc3^3
            tt16 <- 1/xi^4
            tt18 <- dd8 * dd7 - bb1 * aa2 * cc2 * dd6
            l3_10 <- (-jj13 * tt18^3) - 3 * jj13 * (ee1 * 
                aa2 * hh4 * jj08 + 2 * bb1 * dd8 * cc2 * dd6 - 
                2 * jj12 * dd7) * tt18 - jj13 * ((-2 * kk1 * 
                aa2 * qq05 * tt08) - 3 * ee1 * dd8 * hh4 * jj08 - 
                6 * bb1 * jj12 * cc2 * dd6 + 6 * tt16 * dd7) + 
                6 * tt16 * dd3 * dd7 - 6 * jj12 * dd7 - 6 * bb1 * 
                jj12 * dd3 * cc2 * dd6 + 6 * bb1 * dd8 * cc2 * 
                dd6 - 3 * ee1 * dd8 * dd3 * hh4 * jj08 + 3 * 
                ee1 * aa2 * hh4 * jj08 - 2 * kk1 * aa2 * dd3 * 
                qq05 * tt08

        out[["d3"]] <- list(l3_1, l3_2, l3_3, l3_4, l3_5, l3_6, l3_7, l3_8, l3_9, l3_10)

        if (deriv > 3) {
            uu1 <- 1/exp1^(4 * rho)
            uu2 <- xi^3
            l4_1 <- uu1 * (ee3 - 3) * (ee3 - 2) * (ee3 - 1) * 
                uu2 * aa1^(ee3 - 4) + (6 * uu1 * uu2 * (xi + 
                1))/aa1^4
            vv09 <- ee3 - 3
            l4_2 <- 3 * kk1 * ll8 * ff7 * kk2 * ll5^vv09 + 
                uu1 * vv09 * ll8 * ff7 * uu2 * cc2 * ll5^(ee3 - 
                  4) - (6 * kk1 * kk2 * dd3)/ll5^3 + (6 * uu1 * 
                uu2 * dd3 * cc2)/ll5^4
            ww11 <- gg7 - 3
            ww12 <- cc3^(gg7 - 4)
            ww15 <- cc3^ww11
            l4_3 <- (-kk1 * mm11 * oo10 * kk2 * ww15 * (log.cc3/kk2 - 
                (bb1 * aa2 * cc2)/cc3)) + 2 * kk1 * mm11 * xi * 
                ww15 - kk1 * ww15 + uu1 * ww11 * mm11 * kk2 * 
                cc2 * ww12 - uu1 * oo10 * xi * cc2 * ww12 - uu1 * 
                ww11 * xi * cc2 * ww12 + 2 * kk1 * kk2 * tt08 + 
                4 * kk1 * xi * dd3 * tt08 - (6 * uu1 * kk2 * 
                dd3 * cc2)/cc3^4
            l4_4 <- 4 * ee1 * ff7 * xi * ll5^ll8 + 5 * kk1 * 
                ll8 * ff7 * kk2 * cc2 * ll5^vv09 + uu1 * vv09 * 
                ll8 * ff7 * uu2 * hh4 * ll5^(ee3 - 4) + (4 * 
                ee1 * xi * dd3)/ll5^2 - (10 * kk1 * kk2 * dd3 * 
                cc2)/ll5^3 + (6 * uu1 * uu2 * dd3 * hh4)/ll5^4
            yy18 <- log.cc3/kk2
            l4_5 <- (-2 * ee1 * oo10 * xi * mm12 * (bb1 * 
                mm11 * cc2 * dd6 + yy18)) - kk1 * mm11 * oo10 * 
                kk2 * cc2 * ww15 * (bb1 * ww11 * cc2 * dd6 + 
                yy18) - 2 * ee1 * aa2 * mm12 - 2 * ee1 * oo10 * 
                mm12 - 2 * kk1 * mm11 * oo10 * xi * cc2 * ww15 - 
                kk1 * oo10 * cc2 * ww15 - kk1 * mm11 * cc2 * 
                ww15 - 2 * ee1 * dd3 * jj08 - 2 * ee1 * xi * 
                jj08 + 2 * kk1 * kk2 * cc2 * tt08 + 8 * kk1 * 
                xi * dd3 * cc2 * tt08 - (6 * kk2 * dd3 * cc2^2)/(exp1^(4 * 
                rho) * cc3^4)
            l4_6 <- ee1 * oo10 * xi * mm12 * tt18^2 - 2 * 
                ee1 * mm12 * tt18 - 2 * kk1 * mm11 * xi * cc2 * 
                ww15 * tt18 + 2 * kk1 * cc2 * ww15 * tt18 + ee1 * 
                oo10 * xi * mm12 * (ee1 * aa2 * hh4 * jj08 + 
                2 * bb1 * dd8 * cc2 * dd6 - (2 * dd7)/xi^3) + 
                4 * kk1 * cc2 * ww15 + 2 * uu1 * ww11 * xi * 
                hh4 * ww12 - 4 * uu1 * hh4 * ww12 + 2 * ee1 * 
                jj08 - 4 * kk1 * dd3 * cc2 * tt08 - 4 * kk1 * 
                xi * cc2 * tt08 + (6 * uu1 * xi * dd3 * hh4)/cc3^4
            l4_7 <- bb1 * cc3^ff7 + 7 * ee1 * ff7 * xi * cc2 * 
                cc3^ll8 + 6 * kk1 * ll8 * ff7 * kk2 * hh4 * cc3^vv09 + 
                uu1 * vv09 * ll8 * ff7 * uu2 * qq05 * cc3^(ee3 - 
                  4) - (bb1 * dd3)/cc3 + (7 * ee1 * xi * dd3 * 
                cc2)/cc3^2 - (12 * kk1 * kk2 * dd3 * hh4)/cc3^3 + 
                (6 * uu1 * uu2 * dd3 * qq05)/cc3^4
            l4_8 <- (-bb1 * cc3^oo10 * (bb1 * oo10 * cc2 * 
                dd6 + yy18)) - 3 * ee1 * oo10 * xi * cc2 * mm12 * 
                (bb1 * mm11 * cc2 * dd6 + yy18) - kk1 * mm11 * 
                oo10 * kk2 * hh4 * ww15 * (bb1 * ww11 * cc2 * 
                dd6 + yy18) - 3 * ee1 * aa2 * cc2 * mm12 - 3 * 
                ee1 * oo10 * cc2 * mm12 - 2 * kk1 * mm11 * oo10 * 
                xi * hh4 * ww15 - kk1 * oo10 * hh4 * ww15 - kk1 * 
                mm11 * hh4 * ww15 + bb1 * dd6 - 4 * ee1 * dd3 * 
                cc2 * jj08 - 3 * ee1 * xi * cc2 * jj08 + 2 * 
                kk1 * kk2 * hh4 * tt08 + 10 * kk1 * xi * dd3 * 
                hh4 * tt08 - (6 * kk2 * dd3 * cc2^3)/(exp1^(4 * 
                rho) * cc3^4)
            ad17 <- 2 * bb1 * dd8 * cc2 * dd6
            ad19 <- -(2 * dd7)/xi^3
            ad20 <- cc3^oo10
            ad21 <- dd8 * dd7
            ad22 <- ad21 + bb1 * mm11 * cc2 * dd6
            l4_9 <- bb1 * ad20 * (bb1 * oo10 * cc2 * dd6 + 
                ad21)^2 + ee1 * oo10 * xi * cc2 * mm12 * ad22^2 + 
                2 * ee1 * aa2 * cc2 * mm12 * ad22 + 2 * ee1 * 
                oo10 * cc2 * mm12 * ad22 + bb1 * ad20 * ((-ee1 * 
                oo10 * hh4 * jj08) + ad17 + ad19) + ee1 * oo10 * 
                xi * cc2 * mm12 * ((-ee1 * mm11 * hh4 * jj08) + 
                ad17 + ad19) + 4 * ee1 * cc2 * jj08 - 6 * kk1 * 
                dd3 * hh4 * tt08 - 4 * kk1 * xi * hh4 * tt08 + 
                (6 * xi * dd3 * cc2^3)/(exp1^(4 * rho) * cc3^4)
            ae16 <- dd8 * dd7 + bb1 * pp07 * cc2 * dd6
            l4_10 <- (-bb1 * pp08 * ae16^3) - 3 * bb1 * pp08 * 
                ((-ee1 * pp07 * hh4 * jj08) + 2 * bb1 * dd8 * 
                  cc2 * dd6 - 2 * jj12 * dd7) * ae16 - bb1 * 
                pp08 * (2 * kk1 * pp07 * qq05 * tt08 - 3 * ee1 * 
                dd8 * hh4 * jj08 - 6 * bb1 * jj12 * cc2 * dd6 + 
                (6 * dd7)/xi^4) + 6 * kk1 * hh4 * tt08 - (6 * 
                (xi + 1) * qq05)/(exp1^(4 * rho) * cc3^4)
            af05 <- cc2^4
            l4_11 <- bb1 * cc2 * cc3^ff7 + 7 * ee1 * ff7 * 
                xi * hh4 * cc3^ll8 + 6 * kk1 * ll8 * ff7 * kk2 * 
                qq05 * cc3^vv09 + uu1 * vv09 * ll8 * ff7 * uu2 * 
                af05 * cc3^(ee3 - 4) - (bb1 * dd3 * cc2)/cc3 + 
                (7 * ee1 * xi * dd3 * hh4)/cc3^2 - (12 * kk1 * 
                kk2 * dd3 * qq05)/cc3^3 + (6 * uu1 * uu2 * dd3 * 
                af05)/cc3^4
            ag23 <- log.cc3/kk2 - bb1 * aa2 * cc2 * dd6
            l4_12 <- (-bb1 * cc2 * cc3^oo10 * ag23) - 3 * 
                ee1 * oo10 * xi * hh4 * mm12 * ag23 - kk1 * mm11 * 
                oo10 * kk2 * qq05 * ww15 * ag23 + 4 * ee1 * hh4 * 
                mm12 + 5 * kk1 * mm11 * xi * qq05 * ww15 - 4 * 
                kk1 * qq05 * ww15 + uu1 * ww11 * mm11 * kk2 * 
                af05 * ww12 - uu1 * oo10 * xi * af05 * ww12 - 
                uu1 * ww11 * xi * af05 * ww12 + bb1 * cc2 * dd6 - 
                4 * ee1 * dd3 * hh4 * jj08 - 3 * ee1 * xi * hh4 * 
                jj08 + 2 * kk1 * kk2 * qq05 * tt08 + 10 * kk1 * 
                xi * dd3 * qq05 * tt08 - (6 * uu1 * kk2 * dd3 * 
                af05)/cc3^4
            ah24 <- (-(2 * dd7)/xi^3) + 2 * bb1 * dd8 * cc2 * 
                dd6 + ee1 * aa2 * hh4 * jj08
            ah27 <- tt18^2
            l4_13 <- bb1 * cc2 * ad20 * ah27 + ee1 * oo10 * 
                xi * hh4 * mm12 * ah27 - 4 * ee1 * hh4 * mm12 * 
                tt18 - 2 * kk1 * mm11 * xi * qq05 * ww15 * tt18 + 
                2 * kk1 * qq05 * ww15 * tt18 + bb1 * cc2 * ad20 * 
                ah24 + ee1 * oo10 * xi * hh4 * mm12 * ah24 + 
                6 * kk1 * qq05 * ww15 + 2 * uu1 * ww11 * xi * 
                af05 * ww12 - 4 * uu1 * af05 * ww12 + 4 * ee1 * 
                hh4 * jj08 - 6 * kk1 * dd3 * qq05 * tt08 - 4 * 
                kk1 * xi * qq05 * tt08 + (6 * uu1 * xi * dd3 * 
                af05)/cc3^4
            l4_14 <- (-bb1 * cc2 * pp08 * ae16^3) - 3 * bb1 * 
                cc2 * pp08 * ((-ee1 * pp07 * hh4 * jj08) + 2 * 
                bb1 * dd8 * cc2 * dd6 - 2 * jj12 * dd7) * ae16 - 
                bb1 * cc2 * pp08 * (2 * kk1 * pp07 * qq05 * tt08 - 
                  3 * ee1 * dd8 * hh4 * jj08 - 6 * bb1 * jj12 * 
                  cc2 * dd6 + (6 * dd7)/xi^4) + 6 * kk1 * qq05 * 
                tt08 - (6 * (xi + 1) * cc2^4)/(exp1^(4 * rho) * 
                cc3^4)
            aj08 <- 1/cc3^4
            aj20 <- 1/xi^5
            aj23 <- (-2 * jj12 * dd7) + 2 * bb1 * dd8 * cc2 * 
                dd6 + ee1 * aa2 * hh4 * jj08
            l4_15 <- (-jj13 * tt18^4) - 6 * jj13 * aj23 * 
                tt18^2 - 3 * jj13 * aj23^2 - 4 * jj13 * ((-2 * 
                kk1 * aa2 * qq05 * tt08) - 3 * ee1 * dd8 * hh4 * 
                jj08 - 6 * bb1 * jj12 * cc2 * dd6 + 6 * tt16 * 
                dd7) * tt18 - jj13 * (6 * uu1 * aa2 * af05 * 
                aj08 + 8 * kk1 * dd8 * qq05 * tt08 + 12 * ee1 * 
                jj12 * hh4 * jj08 + 24 * bb1 * tt16 * cc2 * dd6 - 
                24 * aj20 * dd7) - 24 * aj20 * dd3 * dd7 + 24 * 
                tt16 * dd7 + 24 * bb1 * tt16 * dd3 * cc2 * dd6 - 
                24 * bb1 * jj12 * cc2 * dd6 + 12 * ee1 * jj12 * 
                dd3 * hh4 * jj08 - 12 * ee1 * dd8 * hh4 * jj08 + 
                8 * kk1 * dd8 * dd3 * qq05 * tt08 - 8 * kk1 * 
                aa2 * qq05 * tt08 + 6 * uu1 * aa2 * dd3 * af05 * 
                aj08

          out[["d4"]] <- list(l4_1, l4_2, l4_3, l4_4, l4_5, l4_6, l4_7, l4_8, l4_9, l4_10,
                              l4_11, l4_12, l4_13, l4_14, l4_15)

        }

      }

    }

  }

  return( out )

}
