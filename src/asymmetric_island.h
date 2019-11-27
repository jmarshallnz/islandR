#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <Rcpp.h>

class logdouble; // forward definition

class Island {
  // Not very efficient 2d,3d array class generated from vectors
  typedef std::vector< std::vector<double> > NumericArray2;
  typedef std::vector< std::vector<std::vector<double> > > NumericArray3;

  int ng;							// # groups
  std::vector< Rcpp::IntegerMatrix > MLST;		// haps for each
  Rcpp::IntegerVector size;				// size of each
  Rcpp::IntegerVector nST;				// # unique STs in each
  NumericArray2 FREQ;		// freq of STs in each group
  NumericArray2 ABUN;	  // abundance of STs in each group
  NumericArray3 acount;	// acount[i][j][k] gives the count, in pop i, and locus j, of allele k

  int nloc;						// # loci
	bool init;

	Rcpp::IntegerMatrix human;				// those sampled from humans

	Rcpp::LogicalMatrix human_unique;
	std::vector< Rcpp::LogicalMatrix > beast_unique;
	bool ****same;
	bool *****ksame;

public:
  // output stuff
  Rcpp::NumericMatrix evolution_traces;
  Rcpp::NumericMatrix accept;
  Rcpp::List human_likelihoods;

  Island() {
		init = false;
		nloc = 7;
		same = NULL;
		ksame = NULL;
	}
  void initialise(Rcpp::IntegerMatrix isolates);

	// mcmc6f infers M and R from seqs of known origin, sampling nsamples after nburnin samples, with given thinning
	void mcmc6f(const double beta, const Rcpp::NumericVector &gamma_m, const Rcpp::NumericVector &gamma_r, const int samples, const int burnin, const int thin);

	~Island() {
		/* free memory */
		for(int i = 0; i < human.nrow(); i++) {
		  for(int ii = 0; ii < ng; ii++) {
		    for(int jj = 0; jj < nST[ii]; jj++) {
		      delete[] same[i][ii][jj];
		    }
		    delete[] same[i][ii];
		  }
		  delete[] same[i];
		}
		delete[] same;

		for(int i = 0; i < ng; i++) {
		  for(int j = 0; j < nST[i]; j++) {
		    for(int ii = 0; ii < ng; ii++) {
		      for(int jj = 0; jj < nST[ii]; jj++) {
		        delete[] ksame[i][j][ii][jj];
		      }
		      delete[] ksame[i][j][ii];
		    }
		    delete[] ksame[i][j];
		  }
		  delete[] ksame[i];
		}
		delete[] ksame;
	}

	int multinom(const Rcpp::NumericVector &p);
	int sample(int n);

	double likHi6(const int id, const int i, const Rcpp::NumericMatrix &A, const NumericArray3 &b, const Rcpp::NumericMatrix &M, const Rcpp::NumericMatrix &R);
	double known_source_loglik(const Rcpp::NumericMatrix &A, const NumericArray3 &b, const Rcpp::NumericMatrix &M, const Rcpp::NumericMatrix &R);
	logdouble known_source_loglik_ij(int i, int j, const Rcpp::NumericMatrix &A, const NumericArray3 &b, const Rcpp::NumericMatrix &M, const Rcpp::NumericMatrix &R);

	NumericArray3 calc_b(const Rcpp::NumericMatrix &A);
	Rcpp::NumericVector normalise(const Rcpp::NumericVector &x);
	Rcpp::NumericMatrix normalise_rows(const Rcpp::NumericMatrix &x);
	Rcpp::NumericMatrix matrix_product(const Rcpp::NumericMatrix &x, const Rcpp::NumericMatrix &y);
	void precalc();
};

#endif//_CLUSTER_H_
