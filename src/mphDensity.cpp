#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////////////////
//  Laplace inversion of MPH density  //
////////////////////////////////////////

//' Inverse Laplace of 2D MPH Laplace transformed density
//'
//' Laplace inverse with Abate-Whitt frame work
//'
//' @param t1 First argument in which to evaluate the density
//' @param t2 Second argument in which to evaluate the density
//' @param eta Constants for Abate-Whitt
//' @param beta Constants for Abate-Whitt
//' @param alpha Inital probability vector
//' @param U Inverse of minus the generator matrix
//' @param R Reward matrix
//'
//' @return MPH density evaluation
// [[Rcpp::export]]

double mphDensity(arma::cx_double t1,
                   arma::cx_double t2,
                   arma::cx_vec eta,
                   arma::cx_vec beta,
                   arma::vec alpha,
                   arma::mat U,
                   arma::mat R){
  int n = alpha.n_rows;
  int m = eta.n_rows;
  arma::cx_double h4 = 0;
  for (unsigned int i = 0; i < m; i++) {
    for (unsigned int j = 0; j < m; j++){
      arma::cx_vec theta = {beta[i]/t1, beta[j]/t2};
      arma::cx_double inverse = arma::accu(alpha.t()*arma::inv(arma::eye(n,n) + U*arma::diagmat(R*theta)));
      h4 = h4 + eta[i]*eta[j]*inverse;
    }
  }
  return real(h4/(t1*t2));
}
