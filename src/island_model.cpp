#include <Rcpp.h>

#include "cluster.h"

using namespace Rcpp;

// run island model

// [[Rcpp::export]]
List island(IntegerMatrix isolates, int niter = 10000, int seed = -5) {
  int run = 0;

  // convert our isolate matrix to the appropriate format
  myutils::Matrix<int> iso(isolates.nrow(), isolates.ncol());

  for (size_t i = 0; i < isolates.nrow(); i++) {
    for (size_t j = 0; j < isolates.ncol(); j++) {
      iso[i][j] = isolates(i,j);
    }
  }

  // create model class
  Cluster clust;
  clust.initialise(iso);

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
  for (std::map<int, myutils::Matrix<double> >::iterator it = clust.human_likelihoods.begin();
       it != clust.human_likelihoods.end(); ++it) {
    NumericMatrix out(it->second.nrows(), it->second.ncols());
    for (int i = 0; i < out.nrow(); i++) {
      for (int j = 0; j < out.ncol(); j++)
        out(i,j) = it->second[i][j];
    }
    std::stringstream s; s << it->first;
    ret[s.str()] = out;
  }
  return ret;
}
