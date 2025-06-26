#ifndef WLS_COV_H
#define WLS_COV_H

#include <RcppArmadillo.h>

// Declaration of the function
Rcpp::List WLS_cov(arma::mat D, arma::mat X, arma::mat W);

#endif
