#include <Rcpp.h>
#include "logdouble.h"

using namespace Rcpp;

// log-likelihood function for assignment model

// [[Rcpp::export]]
double log_lik(NumericMatrix humans, List Phi, NumericVector p) {
  logdouble lik;
  // convert p to exp(p)
  NumericVector log_phi = Phi[0];
  NumericMatrix phi = Phi[1];
  NumericVector P = exp(p);
  P.push_back(1);
  double scaleP = 1 / std::accumulate(P.begin(), P.end(), 0.0);
  for (int j = 0; j < P.size(); j++)
    P[j] = P[j] * scaleP;

  // run through the rows of the human matrix
  for (int h = 0; h < humans.nrow(); ++h) {
    // calculate the likelihood for this human isolate
    double lik_h = 0;
    for (int j = 0; j < P.size(); j++) {
      lik_h += phi(humans(h,0)-1, j) * P[j];
    }
    lik *= logdouble(log_phi[humans(h,0)-1], lik_h) ^ humans(h,1);
  }
  return lik.log();
}
