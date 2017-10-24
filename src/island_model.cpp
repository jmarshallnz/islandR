#include <Rcpp.h>

#include "asymmetric_island.h"

using namespace Rcpp;

// run island model

// [[Rcpp::export]]
List island(IntegerMatrix isolates, NumericVector beta_migration, NumericVector gamma_recombination, int niter = 10000) {

  // create model class
  Island island;
  island.initialise(isolates);

  double beta  = beta_migration[0];
  double gamma = gamma_recombination[0];

  int thin  = 50;

  // output traces. The thin here is how often we sample the evolution parameters.
  // it is NOT how often we sample the probabilities, which is the only thing we care about ATM.
  island.mcmc6f(beta, gamma, niter, thin);

  List out;
  out["hum_lik"] = island.human_likelihoods;
  out["evolution"] = island.evolution_traces;
  return out;
}
