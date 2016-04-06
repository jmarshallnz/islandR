#include <Rcpp.h>

#include "cluster.h"

using namespace Rcpp;

// run island model

// [[Rcpp::export]]
List island(IntegerMatrix isolates, int niter = 10000) {

  // create model class
  Cluster clust;
  clust.initialise(isolates);

  double alpha = 1.0;
  double beta  = 1.0;
  double gamma = 1.0;

  int thin  = 50;

  // output traces. The thin here is how often we sample the evolution parameters.
  // it is NOT how often we sample the probabilities, which is the only thing we care about ATM.
  clust.mcmc6f(alpha, beta, gamma, niter, thin);

  return clust.human_likelihoods;
}
