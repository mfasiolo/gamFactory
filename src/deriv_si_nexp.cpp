// src/deriv_si_nexp.cpp
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// sigmoid and its derivatives
inline double sigmoid(double u){
  if(u >= 0){
    double e = std::exp(-u);
    return 1.0/(1.0+e);
  } else {
    double e = std::exp(u);
    return e/(1.0+e);
  }
}

// Flatten upper triangle (j<=k) into column index (0-based)
inline int uptri_col(int j, int k, int m){
  return j*m - j*(j-1)/2 + (k - j);
}

// [[Rcpp::export(name = ".deriv_si_nexp_cpp")]]
Rcpp::List deriv_si_nexp_cpp(const NumericMatrix& X_si,
                             const NumericMatrix& X_nexp,
                             const NumericVector& param,
                             int deriv,
                             Rcpp::Nullable<NumericVector> alpha_center,
                             double Z0,
                             bool positive_si = false) { 
  
  if(deriv > 3){ Rcpp::warning("deriv exceeds 3; truncating to 3"); deriv = 3; }
  
  const int n      = X_si.nrow();
  const int n_si   = X_si.ncol();     // p
  const int n_nexp = X_nexp.ncol();   // q
  const int m_par  = n_si + n_nexp;   // m
  
  if(param.size() != m_par){
    stop("param length must be q + p (alpha_nexp then alpha_si).");
  }
  
  // Split parameters: param = [alpha_nexp (q), alpha_si (p)]
  NumericVector a_ne(n_nexp), a_si(n_si);
  for(int j=0;j<n_nexp;++j) a_ne[j] = param[j];
  for(int j=0;j<n_si;++j) {
    if (positive_si) {
      a_si[j] = std::exp(param[n_nexp + j]);
    } else {
      a_si[j] = param[n_nexp + j];
    }
  }
  
  // Optional alpha_center
  bool use_center = alpha_center.isNotNull();
  NumericVector a_center;
  if(use_center){
    a_center = alpha_center.get();
    if(a_center.size() != n_si) stop("alpha_center must have length = ncol(X_si).");
  }
  
  // Compute z_t, u_t, om_t first
  NumericVector z(n), u(n), om(n);
  for(int t=0;t<n;++t){
    double zt = 0.0;
    for(int j=0;j<n_si;++j){
      double coef = a_si[j] + (use_center ? a_center[j] : 0.0);
      zt += X_si(t,j) * coef;
    }
    z[t] = zt;
    
    double ut = 0.0;
    for(int j=0;j<n_nexp;++j) ut += X_nexp(t,j) * a_ne[j];
    u[t]  = ut;
    om[t] = sigmoid(ut);
  }
  
  // Output containers
  NumericVector d0(n);                       // s_t
  NumericMatrix d1;                          // n x (q+p)
  NumericMatrix d2;                          // n x choose(m_par+1,2)
  NumericMatrix d3;                          // n x choose(m_par+2,3)
  
  const bool need1 = (deriv >= 1);
  const bool need2 = (deriv >= 2);
  const bool need3 = (deriv >= 3);
  
  if(need1) d1 = NumericMatrix(n, m_par);
  if(need2) d2 = NumericMatrix(n, m_par*(m_par+1)/2);
  if(need3) d3 = NumericMatrix(n, (m_par*(m_par+1)*(m_par+2))/6);
  
  // Variables from previous time step
  double s_prev;
  
  // 1st-order: e = ds/d alpha_si (p), d = ds/d alpha_nexp (q)
  NumericVector e_prev(n_si);      // t-1
  NumericVector d_prev(n_nexp);    // t-1
  
  // 2nd-order
  NumericMatrix h_prev(n_si, n_nexp);  
  NumericMatrix g_prev(n_nexp, n_nexp);
  
  // 3rd-order
  std::vector<double> T3_prev, T3_curr, U3_prev, U3_curr;
  if(need3){
    T3_prev.assign(n_si * n_nexp * n_nexp, 0.0);
    T3_curr.assign(n_si * n_nexp * n_nexp, 0.0);
    U3_prev.assign(n_nexp * n_nexp * n_nexp, 0.0);
    U3_curr.assign(n_nexp * n_nexp * n_nexp, 0.0);
  }
  
  // --- Initialization for t = 0 ---
  if (R_IsNA(Z0)) {
    d0[0] = z[0];
  } else {
    d0[0] = om[0] * Z0 + (1.0 - om[0]) * z[0];
  }
  s_prev = d0[0];
  
  // initialize first derivative and get final_e
  NumericVector e_final_0(n_si); 
  for(int j=0;j<n_si;++j) {
    e_prev[j] = X_si(0,j);
    e_final_0[j] = positive_si ? (e_prev[j] * a_si[j]) : e_prev[j];
  }
  for(int j=0;j<n_nexp;++j) d_prev[j] = 0.0;
  
  if(need1){
    for(int j=0;j<n_nexp;++j) d1(0, j) = 0.0;
    for(int j=0;j<n_si;++j)   d1(0, n_nexp + j) = e_final_0[j];
  }
  
  // Filling in the diagonal elements of the second and third derivatives at t=0
  if (need2 && positive_si) {
    int col = 0;
    for (int j = 0; j < m_par; ++j) {
      for (int k = j; k < m_par; ++k) {
        if (j >= n_nexp && j == k) d2(0, col) = e_final_0[j - n_nexp];
        col++;
      }
    }
  }
  if (need3 && positive_si) {
    int col3 = 0;
    for (int j = 0; j < m_par; ++j) {
      for (int k = j; k < m_par; ++k) {
        for (int l = k; l < m_par; ++l) {
          if (j >= n_nexp && j == k && k == l) d3(0, col3) = e_final_0[j - n_nexp];
          col3++;
        }
      }
    }
  }
  
  // --- Main loop t = 1..n-1 ---
  for (int t = 1; t < n; ++t) {
    const double w  = om[t];
    const double s1 = w * (1.0 - w);                 
    const double s2 = s1 * (1.0 - 2.0 * w);          
    const double s3 = s1 * (1.0 - 6.0*w + 6.0*w*w);  
    
    const double zt = z[t];
    const double s_t = w * s_prev + (1.0 - w) * zt;
    d0[t] = s_t;
    
    // ========= 1st-order =========
    NumericVector e_curr(n_si);
    NumericVector e_final(n_si); // first derivative after mapping
    for (int j = 0; j < n_si; ++j) {
      e_curr[j] = w * e_prev[j] + (1.0 - w) * X_si(t, j);
      e_final[j] = positive_si ? (e_curr[j] * a_si[j]) : e_curr[j]; // chain rule
    }
    
    NumericVector d_curr(n_nexp);
    for (int j = 0; j < n_nexp; ++j) {
      d_curr[j] = w * d_prev[j] + (s_prev - zt) * s1 * X_nexp(t, j);
    }
    
    if (need1) {
      for (int j = 0; j < n_nexp; ++j) d1(t, j) = d_curr[j];
      for (int j = 0; j < n_si;   ++j) d1(t, n_nexp + j) = e_final[j]; // insert mapped e
    }
    
    // ========= 2nd-order =========
    NumericMatrix h_curr(n_si,   n_nexp);
    NumericMatrix g_curr(n_nexp, n_nexp);
    
    if (need2 || need3) {
      for (int i = 0; i < n_si; ++i) {
        const double term_i = e_prev[i] - X_si(t, i);
        for (int a = 0; a < n_nexp; ++a) {
          h_curr(i, a) = w * h_prev(i, a) + term_i * (s1 * X_nexp(t, a));
        }
      }
      for (int a = 0; a < n_nexp; ++a) {
        for (int b = 0; b < n_nexp; ++b) {
          g_curr(a, b) = w * g_prev(a, b)
          + s1 * ( X_nexp(t, a) * d_prev[b] + d_prev[a] * X_nexp(t, b) )
          + (s_prev - zt) * s2 * (X_nexp(t, a) * X_nexp(t, b));
        }
      }
      
      if (need2) {
        int col = 0;
        // 1. g
        for (int j = 0; j < n_nexp; ++j) {
          for (int k = j; k < n_nexp; ++k) d2(t, col++) = g_curr(j, k);
          // 2. h-matrix cross-block (t(h))
          for (int k = n_nexp; k < m_par; ++k) {
            int i = k - n_nexp;
            double h_val = h_curr(i, j);
            if (positive_si) h_val *= a_si[i]; // Chain Rule: Mixed Second-Order Derivation Scaling
            d2(t, col++) = h_val; 
          }
        }
        // 3. The second-order derivative matrix block H_sisi for pure si parameters
        for (int j = n_nexp; j < m_par; ++j) {
          for (int k = j; k < m_par; ++k) {
            double val = 0.0;
            // <-- Chain rule: Substitute the second-order diagonal of the derivative of exp 
            // (i.e. the e_final we calculated earlier)
            if (positive_si && j == k) val = e_final[j - n_nexp];
            d2(t, col++) = val;
          }
        }
      }
    }
    
    // ========= 3rd-order =========
    if (need3) {
      std::fill(T3_curr.begin(), T3_curr.end(), 0.0);
      std::fill(U3_curr.begin(), U3_curr.end(), 0.0);
      
      for (int i = 0; i < n_si; ++i) {
        const double term_i = e_prev[i] - X_si(t, i);
        for (int a = 0; a < n_nexp; ++a) {
          for (int b = 0; b < n_nexp; ++b) {
            const int idx = (i*n_nexp + a)*n_nexp + b;
            T3_curr[idx] = w * T3_prev[idx]
            + s1 * ( X_nexp(t, a) * h_prev(i, b) + X_nexp(t, b) * h_prev(i, a) )
              + term_i * (s2 * X_nexp(t, a) * X_nexp(t, b));
          }
        }
      }
      
      for (int c = 0; c < n_nexp; ++c) {
        for (int a = 0; a < n_nexp; ++a) {
          for (int b = 0; b < n_nexp; ++b) {
            const int idx = (a*n_nexp + b)*n_nexp + c;
            double val = w * U3_prev[idx];
            val += (s1 * X_nexp(t, c)) * g_prev(a, b);
            val += s2 * X_nexp(t, c) * ( X_nexp(t, a) * d_prev[b] + d_prev[a] * X_nexp(t, b) );
            val += s1 * ( X_nexp(t, a) * g_prev(b, c) + g_prev(a, c) * X_nexp(t, b) );
            val += d_prev[c] * s2 * (X_nexp(t, a) * X_nexp(t, b));
            val += (s_prev - zt) * s3 * X_nexp(t, c) * (X_nexp(t, a) * X_nexp(t, b));
            U3_curr[idx] = val;
          }
        }
      }
      
      // Assemble d3[t, ]
      int col3 = 0;
      for (int j = 0; j < m_par; ++j) {
        for (int k = j; k < m_par; ++k) {
          for (int l = k; l < m_par; ++l) {
            const int si_cnt = (j >= n_nexp) + (k >= n_nexp) + (l >= n_nexp);
            double val = 0.0;
            
            if (si_cnt == 0) {
              val = U3_curr[(j * n_nexp + k) * n_nexp + l];
            } 
            else if (si_cnt == 1) {
              int idxs[3] = {j, k, l};
              int si_pos = -1; 
              for(int t2=0; t2<3; ++t2) if(idxs[t2] >= n_nexp) { si_pos=t2; break; }
              const int si_idx = idxs[si_pos] - n_nexp; 
              int a_idx = -1, b_idx = -1, cnt = 0;
              for(int t2=0; t2<3; ++t2){
                if(t2 == si_pos) continue;
                if(cnt == 0){ a_idx = idxs[t2]; cnt = 1; } else { b_idx = idxs[t2]; }
              }
              if(a_idx > b_idx) std::swap(a_idx, b_idx);
              val = T3_curr[(si_idx * n_nexp + a_idx) * n_nexp + b_idx];
              
              // Scaling only occurs if `positive_si` is enabled; otherwise, the original T3 is retained.
              if (positive_si) val *= a_si[si_idx]; 
            } 
            else if (si_cnt == 2) {
              // On a linear scale, the second derivative of `alpha_si` is zero, so it is non-zero only when `positive_si`
              if (positive_si) {
                int idxs[3] = {j, k, l};
                int si_idx1 = -1, si_idx2 = -1, a_idx = -1;
                int c_si = 0;
                for(int t2=0; t2<3; ++t2){
                  if(idxs[t2] >= n_nexp) {
                    if(c_si==0) { si_idx1 = idxs[t2] - n_nexp; c_si++; }
                    else { si_idx2 = idxs[t2] - n_nexp; }
                  } else { a_idx = idxs[t2]; }
                }
                if (si_idx1 == si_idx2) val = h_curr(si_idx1, a_idx) * a_si[si_idx1];
              }
            } 
            else if (si_cnt == 3) {
              // In the linear regime, the third derivative of α_si is zero
              if (positive_si) {
                int si_idx1 = j - n_nexp, si_idx2 = k - n_nexp, si_idx3 = l - n_nexp;
                if (si_idx1 == si_idx2 && si_idx2 == si_idx3) val = e_final[si_idx1];
              }
            }
            d3(t, col3++) = val;
          }
        }
      }
      T3_prev.swap(T3_curr);
      U3_prev.swap(U3_curr);
    }
    
    // Update previous step
    if (need2 || need3) {
      h_prev = h_curr;
      g_prev = g_curr;
    }
    e_prev = e_curr;
    d_prev = d_curr;
    s_prev = s_t;
  }
  
  Rcpp::List out = Rcpp::List::create(_["d0"] = d0);
  if(need1) out["d1"] = d1;
  if(need2) out["d2"] = d2;
  if(need3) out["d3"] = d3;
  return out;
}