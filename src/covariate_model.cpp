#include <Rcpp.h>
#include "logdouble.h"

using namespace Rcpp;

// log-likelihood function for assignment model

// [[Rcpp::export]]
double log_lik(NumericMatrix humans, NumericMatrix phi, NumericVector p) {
  logdouble lik;
  // convert p to exp(p)
  NumericVector P = exp(p);
  P.push_back(1);
  double scaleP = 1 / std::accumulate(P.begin(), P.end(), 0.0);
  for (int j = 0; j < P.size(); j++)
    P[j] = P[j] * scaleP;

  // run through the rows of the human matrix
  for (int h = 0; h < humans.nrow(); ++h) {
    // calculate the likelihood for this human isolate
    std::vector<logdouble> lik_h(P.size());
    for (int j = 0; j < P.size(); j++) {
      lik_h[j] = logdouble(phi(humans(h,0)-1, j), 1) * P[j];
    }
    lik *= sum(lik_h) ^ humans(h,1);
  }
  return lik.log();
}
