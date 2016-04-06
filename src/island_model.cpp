#include <Rcpp.h>

#include "cluster.h"

using namespace Rcpp;

// run island model

// [[Rcpp::export]]
List island(IntegerMatrix isolates, int niter = 10000, int seed = -5) {

  // create model class
  Cluster clust;
  clust.initialise(isolates);

  double alpha = 1.0;
  double beta  = 1.0;
  double gamma = 1.0;

  int thin  = 50;

  myutils::Random ran;
  ran.setseed(seed);

  // output traces. The thin here is how often we sample the evolution parameters.
  // it is NOT how often we sample the probabilities, which is the only thing we care about ATM.
  clust.mcmc6f(alpha, beta, gamma, niter, thin, ran);

  // convert our map of output shit into a list
  List ret;
  for (std::map<int, NumericMatrix>::iterator it = clust.human_likelihoods.begin();
       it != clust.human_likelihoods.end(); ++it) {
    std::stringstream s; s << it->first;
    ret[s.str()] = it->second;
  }
  return ret;
}
