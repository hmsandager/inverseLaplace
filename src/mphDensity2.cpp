// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
// [[Rcpp::export]]

double mphDensity2(cx_double t1,
               cx_double t2, 
               cx_vec eta, 
               cx_vec beta, 
               vec alpha, 
               mat U, 
               mat R){
  int n = alpha.n_rows;
  int m = eta.n_rows;
  cx_double h4 = 0;
  for (unsigned int i = 0; i < m; i++) {
    for (unsigned int j = 0; j < m; j++){
      cx_vec theta = {beta[i]/t1, beta[j]/t2};
      cx_double inverse = accu(alpha.t()*inv(eye(n,n) + U*diagmat(R*theta)));
      h4 = h4 + eta[i]*eta[j]*inverse;
    }
  }
  return real(h4/(t1*t2));
}
