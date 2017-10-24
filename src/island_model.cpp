#include <Rcpp.h>

#include "asymmetric_island.h"

using namespace Rcpp;

// run island model

// [[Rcpp::export]]
List island(IntegerMatrix isolates, int niter = 10000) {

  // create model class
  Island island;
  island.initialise(isolates);

  double alpha = 1.0;
  double beta  = 1.0;
  double gamma = 1.0;

  int thin  = 50;

  // output traces. The thin here is how often we sample the evolution parameters.
  // it is NOT how often we sample the probabilities, which is the only thing we care about ATM.
  island.mcmc6f(beta, gamma, niter, thin);

  return island.human_likelihoods;
}
