#include <Rcpp.h>

#include "asymmetric_island.h"

using namespace Rcpp;

// run island model

// [[Rcpp::export]]
List island(IntegerMatrix isolates, NumericVector beta_migration, NumericVector gamma_mutation, NumericVector gamma_recombination, int samples = 100, int burnin = 10, int thin = 100) {

  // create model class
  Island island;
  island.initialise(isolates);

  double beta  = beta_migration[0];

  island.mcmc6f(beta, gamma_mutation, gamma_recombination, samples, burnin, thin);

  List out;
  out["hum_lik"] = island.human_likelihoods;
  out["evolution"] = island.evolution_traces;
  out["accept"] = island.accept;
  return out;
}
