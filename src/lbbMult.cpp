#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat lbbMult(arma::mat x, 
                  arma::vec lpstart,
                  arma::vec lpend, 
                  arma::mat lnn,
                  arma::imat i2) {
  
  arma::mat lbb(x.n_cols, x.n_cols);
  int m = lpend.size() - 1;
  
  for (int rr = 0; rr <= m; rr++) {
    for (int ss = rr; ss <= m; ss++) {
      arma::vec lnn_rs = lnn.col(i2(rr, ss) - 1);
      lbb(arma::span(lpstart[rr], lpend[rr]),
          arma::span(lpstart[ss], lpend[ss])) = 
            x(arma::span::all, 
              arma::span(lpstart[rr], lpend[rr])).t() * 
                (x(arma::span::all, 
                   arma::span(lpstart[ss], lpend[ss])).each_col() % 
                lnn_rs);
      
      if (ss > rr) {
        lbb(arma::span(lpstart[ss], lpend[ss]),
            arma::span(lpstart[rr], lpend[rr])) = 
              lbb(arma::span(lpstart[rr], lpend[rr]),
                  arma::span(lpstart[ss], lpend[ss])).t();
      }
    }
  }
  
  return lbb;
  
}

