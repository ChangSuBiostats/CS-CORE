#ifndef WLS_MEAN_H
#define WLS_MEAN_H

#include <RcppArmadillo.h>

// Declaration of the function
arma::mat WLS_mean(arma::mat D, arma::mat X, arma::mat W);

#endif
