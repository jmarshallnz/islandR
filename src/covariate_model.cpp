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
      // we could speed this up by doing the magic offset stuff first in phi,
      // i.e. storing it such that the maximal exponent is factored out already.
      // The P's are unlikely to be so small as to affect this very much then,
      // so we could do it all with doubles, then log at the end
      lik_h[j] = logdouble(phi(humans(h,0)-1, j), P[j]);
    }
    lik *= sum(lik_h) ^ humans(h,1);
  }
  return lik.log();
}
