#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <Rcpp.h>
#include "myutils.h"
#include <map>

class Cluster {
	int ng;							// # groups
  myutils::Vector< myutils::Matrix<int> > MLST;		// haps for each
  myutils::Vector<int> size;				// size of each
  myutils::Vector<int> nST;				// # unique STs in each
  myutils::Matrix<int> nalleles;			// nalleles[i][j] # unique alleles in group i at locus j
  myutils::Vector< myutils::Vector<double> > FREQ;		// freq of STs in each group
  myutils::Vector< myutils::Vector<double> > ABUN;		// abundance of STs in each group
  myutils::Matrix< myutils::Vector<double> > acount;	// acount[i][j][k] gives the count, in pop i, and locus j, of allele k

  int nloc;						// # loci
	bool init;

	myutils::Matrix<int> human;				// those sampled from humans

	myutils::Matrix<bool> human_unique;
	myutils::Vector< myutils::Matrix<bool> > beast_unique;
	bool ****same;
	bool *****ksame;

public:
  // output stuff
  myutils::Matrix<double> evolution_traces;
	std::map<int, myutils::Matrix<double> > human_likelihoods;

	Cluster() {
		init = false;
		nloc = 7;
		same = NULL;
		ksame = NULL;
	}
  void initialise(myutils::Matrix<int> &isolates);

	// mcmc6f infers M and R from seqs of known origin, and runs 100 side-chains to infer F given M and R
	void mcmc6f(const double alpha, const double beta, const double gamma_, const int niter, const int thin, myutils::Random &ran);

	~Cluster() {
		if(init) {
			int i,j;
			for(i=0;i<ng;i++) {
				for(j=0;j<nloc;j++) {
					acount[i][j].resize(0);
				}
				MLST[i].resize(0,0);
			}
		}
		/* free memory */
		for(int i = 0; i < human.nrows(); i++) {
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

	int multinom(myutils::Vector<double> &p, myutils::Random &ran);
	double likHi6(const int id, const int i, const myutils::Matrix<double> &a, const myutils::Matrix< myutils::Vector<double> > &b, const myutils::Matrix<double> &r);
	double known_source_loglik(const myutils::Matrix<double> &a, const myutils::Matrix< myutils::Vector<double> > &b, const myutils::Matrix<double> &r);

	void recalc_b(myutils::Matrix<double> &a, myutils::Matrix< myutils::Vector<double> > &b);
	void calc_A(Rcpp::NumericMatrix &a, myutils::Matrix<double> &A);
	void calc_Ai(const Rcpp::NumericVector &a, myutils::Matrix<double> &A, const int i);
	void calc_R(myutils::Matrix<double> &r, myutils::Matrix<double> &R);
	void calc_Ri(myutils::Matrix<double> &r, myutils::Matrix<double> &R, const int i);
	void precalc();
	template<typename T> void pSWAP(myutils::Vector<T> &a, myutils::Vector<T> &b);
	template<typename T> void pSWAP(myutils::Matrix<T> &a, myutils::Matrix<T> &b);
};

#endif//_CLUSTER_H_
