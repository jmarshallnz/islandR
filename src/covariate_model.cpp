#include <Rcpp.h>
using namespace Rcpp;

// log-likelihood function for assignment model

// [[Rcpp::export]]
double log_lik(NumericMatrix humans, NumericMatrix phi, NumericVector p) {
  double loglik = 0;

  // run through the rows of the human matrix
  for (size_t h = 0; h < humans.nrow(); ++h) {
    // calculate the likelihood for this human isolate
    double lik_h = 0;
    for (size_t j = 0; j < p.size(); j++) {
      lik_h += p[j] * phi(humans(h,0)-1, j);
    }
    loglik += humans(h,1) * log(lik_h);
  }
  return loglik;
}
