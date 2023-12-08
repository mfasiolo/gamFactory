#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

colvec invlogit(colvec x) {
  int n = x.size();
  colvec result(n);
  
  for (int i=0; i < n; ++i) {
    result[i] = 1.0 / (1.0 + exp (-1.0 * x[i]));
  }
  return result;
}

// [[Rcpp::export(.expsmooth_cpp)]]
List expsmooth(NumericVector y, NumericMatrix Xi, NumericVector beta, 
                   NumericVector x0, int deriv) {
  
  int n = y.length();
  int p = beta.length();
  
  mat XiMat = as<mat>(Xi);
  colvec betaCol = as<colvec>(beta);
  
  colvec w = invlogit(XiMat * betaCol);
  
  // Perform smoothing
  colvec g(n);
  g[0] = as<double>(x0) * w[0] + (1 - w[0]) * y[0];
  
  for (int ii = 1; ii < n; ++ii) {
    g[ii] = w[ii] * g[ii - 1] + (1 - w[ii]) * y[ii];
  }
  
  // Derivatives of g w.r.t. beta
  mat d1, d2, d3;
  
  if (deriv) {
    colvec a1 = w % (1 - w);
    mat w1 = XiMat;
    w1.each_col() %= a1;
    
    colvec gmx(n);
    gmx(0) = x0[0];
    gmx(span(1,n-1)) = g(span(0,n-2));
    gmx -= y; 
    
    d1 = mat(n, p);
    
    d1.row(0) = w1.row(0) * gmx[0];
    
    for (int ii = 1; ii < n; ++ii) {
      d1.row(ii) = w1.row(ii) * gmx[ii] + w[ii] * d1.row(ii - 1);
    }
    
    if (deriv > 1) {
      colvec a2 = a1 - 2 * a1 % w;
      mat w2 = mat(n, p * (p + 1) / 2);
      d2 = mat(n, p * (p + 1) / 2);
      
      int zz = 0;
      for (int jj = 0; jj < p; ++jj) {
        for (int kk = jj; kk < p; ++kk) {
          w2.col(zz) = a2 % XiMat.col(jj) % XiMat.col(kk);
          d2.row(0).col(zz) = w2(0, zz) * gmx[0];
          
          for (int ii = 1; ii < n; ++ii) {
            d2.row(ii).col(zz) = w2(ii, zz) * gmx[ii] + w1(ii, jj) * d1.row(ii - 1).col(kk) +
              w1(ii, kk) * d1.row(ii - 1).col(jj) + w[ii] * d2.row(ii - 1).col(zz);
          }
          
          ++zz;
        }
      }
      
      if (deriv > 2) {
        colvec a3 = a2 - 2 * a2 % w - 2 * square(a1);
        mat w3 = mat(n, p * (p + 1) * (p + 2) / 6);
        d3 = mat(n, p * (p + 1) * (p + 2) / 6);
        
        int zz = 0;
        mat ind(p, p);
        for (int ir = 0; ir < p; ++ir) {
          for (int ic = ir; ic < p; ++ic) {
            ind(ir, ic) = zz;
            ++zz;
          }
        }
        
        zz = 0;
        for (int jj = 0; jj < p; ++jj) {
          for (int kk = jj; kk < p; ++kk) {
            for (int ll = kk; ll < p; ++ll) {
              w3.col(zz) = a3 % XiMat.col(jj) % XiMat.col(kk) % XiMat.col(ll);
              d3.row(0).col(zz) = w3(0, zz) * gmx[0];
              
              for (int ii = 1; ii < n; ++ii) {
                d3.row(ii).col(zz) = w3(ii, zz) * gmx[ii] +
                  w2(ii, ind(jj, kk)) * d1.row(ii - 1).col(ll) +
                  w2(ii, ind(jj, ll)) * d1.row(ii - 1).col(kk) +
                  w2(ii, ind(kk, ll)) * d1.row(ii - 1).col(jj) +
                  w1(ii, jj) * d2.row(ii - 1).col(ind(kk, ll)) +
                  w1(ii, kk) * d2.row(ii - 1).col(ind(jj, ll)) +
                  w1(ii, ll) * d2.row(ii - 1).col(ind(jj, kk)) +
                  w[ii] * d3.row(ii - 1).col(zz);
              }
              
              ++zz;
            }
          }
        }
      }
    }
  }
  
  return List::create(Named("d0", g), Named("d1", d1), Named("d2", d2), Named("d3", d3));
}