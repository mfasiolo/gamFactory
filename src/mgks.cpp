#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.mgks_cpp)]]
List mgks(NumericMatrix y, arma::field<arma::mat> dist, NumericVector beta, int deriv) {
  
  mat yA = as<mat>(y);
  colvec betaA = as<colvec>(beta);
  rowvec w = exp(2 * beta);
  
  int n = dist(0).n_rows;
  int n0 = dist(0).n_cols;
  int p = beta.length();
  
  std::string case_type;
  colvec d0(n);
  mat d1, d2, d3;
  mat DaDb, D2alD2b, ind2, DaDbKK_DlkDb_c, D3alD3b;
  colvec DaDb_DlkDb;
  if (deriv) {
    d1 = mat(n, p);
    if (deriv > 1) {
      DaDb = mat(n0, p);
      D2alD2b = mat(n0, p * (p + 1) / 2);
      DaDb_DlkDb = colvec(p * (p + 1) / 2);
      d2 = mat(n, p * (p + 1) / 2);
      if (deriv > 2) {
        ind2 = mat(p, p);
        d3 = mat(n, Rf_choose(p + 2, p - 1));
        int zz = 0;
        for (int ir = 0; ir < p; ++ir) {
          for (int ic = ir; ic < p; ++ic) {
            ind2(ir, ic) = zz;
            zz++;
          }
        }
        DaDbKK_DlkDb_c = mat(n0, Rf_choose(p + 2, p - 1));
        D3alD3b = mat(n0, Rf_choose(p + 2, p - 1));
      }
    }
  }
  
  colvec yii = yA.col(0);
  mat dist_ii(n0, p);
  for (int ii = 0; ii < n; ++ii) {
    
    if (yA.n_cols > 1) {
      yii = yA.row(ii).t();
    } 
    
    for(int kk = 0; kk < p; ++kk){
     dist_ii.col(kk) = dist(kk).row(ii).t();
    }
    
    // dist[j, k] = (X0[j, k] - xi[k])^2 * w[k] for j = 1, ..., n0 and k = 1, ..., d
    // dist <- - t((tX0 - xi)^2 * w)
    dist_ii.each_row() %= -w;
    
    // Vector of log-kernels logK[j] = sum_k (X0[j, k] - xi[k])^2 * w[k] for j = 1, ..., n0 
    colvec logK = sum(dist_ii, 1);
    
    // sum exp trick
    double mx = max(logK);
    colvec al = exp(logK - mx) / sum(exp(logK - mx));
    
    d0(ii) = sum(al % yii);
    
    // Derivatives of g w.r.t. beta
    if (deriv) {
      mat DlkDb = 2 * dist_ii;  // 2 * dist[kk, jj] = D logK[kk] / D beta[jj]  
      colvec sum_al_DlkDb = sum(DlkDb.each_col() % al, 0).t();
      mat DlkDb_c = DlkDb.each_row() - sum_al_DlkDb.t();
      
      
      // DaDb[i,j] = Dalpha[i] / D beta[j]
      for (int jj = 0; jj < p; ++jj) {
        colvec DaDbJJ = al % DlkDb_c.col(jj);
        d1(ii, jj) = sum(yii % DaDbJJ);
        // Store stuff needed for higher derivs
        if (deriv > 1) {
          DaDb.col(jj) = DaDbJJ;
        }
      }
      
      if (deriv > 1) {
        int zz = 0;
        for (int jj = 0; jj < p; ++jj) {
          for (int kk = jj; kk < p; ++kk) {
            DaDb_DlkDb(zz) = sum(DaDb.col(kk) % DlkDb.col(jj));
            colvec DaDbKK_DlkDb_cJJ = DaDb.col(kk) % DlkDb_c.col(jj);
            D2alD2b.col(zz) = DaDbKK_DlkDb_cJJ - al * DaDb_DlkDb(zz);
            if (jj == kk) {  // Diagonal entries need an extra term
              D2alD2b.col(zz) += al % (2 * DlkDb_c.col(jj));
            }
            // Store stuff needed for higher derivs
            if (deriv > 2) {
              DaDbKK_DlkDb_c.col(zz) = DaDbKK_DlkDb_cJJ; 
            }
            d2(ii, zz) = sum(yii % D2alD2b.col(zz));
            zz++;
          }
        }
        
        if (deriv > 2) {
          zz = 0;
          for (int jj = 0; jj < p; ++jj) {
            for (int kk = jj; kk < p; ++kk) {
              for (int ll = kk; ll < p; ++ll) {
                D3alD3b.col(zz) = D2alD2b.col(ind2(kk, ll)) % DlkDb_c.col(jj) -
                  DaDb.col(kk) * DaDb_DlkDb(ind2(jj, ll)) -
                  DaDb.col(ll) * DaDb_DlkDb(ind2(jj, kk)) -
                  al * sum(D2alD2b.col(ind2(kk, ll)) % DlkDb.col(jj));
                
                if (jj == kk) {
                  D3alD3b.col(zz) += 2 * DaDbKK_DlkDb_c.col(ind2(jj, ll)) -
                    al * (2 * DaDb_DlkDb(ind2(jj, ll)));
                }
                if (jj == ll) {
                  D3alD3b.col(zz) += 2 * DaDbKK_DlkDb_c.col(ind2(jj, kk)) -
                    al * (2 * DaDb_DlkDb(ind2(jj, kk)));
                  if (jj == kk) {
                    D3alD3b.col(zz) += al % (4 * DlkDb_c.col(jj));
                  }
                }
                d3(ii, zz) = sum(yii % D3alD3b.col(zz));
                zz++; 
              }
            }
          }
        } 
      } 
    }
  }
  
  return List::create(Named("d0", d0), Named("d1", d1), Named("d2", d2), Named("d3", d3));
}