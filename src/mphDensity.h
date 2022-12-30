#include <RcppArmadillo.h>

double mphDensity(arma::cx_double t1,
                   arma::cx_double t2,
                   arma::cx_vec eta,
                   arma::cx_vec beta,
                   arma::vec alpha,
                   arma::mat U,
                   arma::mat R);
