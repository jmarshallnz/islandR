#include "cluster.h"

#include <sstream> // for the stream stuff
#include <iostream>

#include <Rcpp.h>

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::_;
using Rcpp::clone;

using namespace myutils;

double Cluster::known_source_loglik(const NumericMatrix &A, const Array3 &b, const NumericMatrix &R) {
	double loglik = 0.0;
	/* Cycle through each unique ST in each group, taking account of abundance of the STs */
	for (int i = 0; i < ng; i++) {
		double punique = A(i,ng);
		for (int j = 0; j < nST[i]; j++) {
			double ncopiesj = ABUN[i][j];
			std::vector<double> psame(nloc);
			std::vector<double> pdiff(nloc);
			for (int l = 0; l < nloc; l++) {
				int allele = MLST[i][j][l];
				double ac = acount[i][l][allele];
				double ac_ = (ac*(double)size[i]-1.0)/(double)(size[i]-1);
				double bk = b[i][l][allele] - A(i,i)*ac + A(i,i)*ac_;
				double b_ = R(i,0) * bk + R(i,1) * (1.0-A(i,ng));// if no rec then must be same as CF (barring mutation)
				if(fabs(b_)<1.0e-7) {
					b_ = 0.0;
				}
				psame[l] = b_;
				b_ = R(i,0) * bk;	// different so must have been recombination

				if(fabs(b_)<1.0e-7) {
					b_ = 0.0;
				}
				pdiff[l] = b_;
			}
			double l_j = 0.0;
			for (int ii = 0; ii < ng; ii++) {						//	Cycle through source of the clonal frame
				double l_ii = 0.0;
				double mii = A(i,ii)/(1.0-A(i,ng));
				for (int jj = 0; jj < nST[ii]; jj++) {				//	Cycle through each ST from that source
					double ncopiesjj = (i==ii && j==jj) ? ABUN[ii][jj]-MIN(ABUN[ii][jj],1.0)
						: ABUN[ii][jj];
					double l_jj = mii;
					bool *BEAST_UNIQUE = beast_unique[i][j];
					bool *SAME = ksame[i][j][ii][jj];
					for (int l = 0; l < nloc; l++, BEAST_UNIQUE++, SAME++) {
						if(*BEAST_UNIQUE) {				// new allele (allow some rounding error)
							l_jj *= punique;
						}
						else if(*SAME) {				// previously observed and same as CF
							l_jj *= psame[l];
						}
						else {							// previously observed but different to CF
							l_jj *= pdiff[l];
						}
					}
					l_ii += l_jj * ncopiesjj;
				}
				l_j += l_ii / size[ii];
			}
			loglik += log(l_j) * ncopiesj;
		}
	}
	return loglik;
}

double Cluster::likHi6(const int id, const int i, const NumericMatrix &A, const Array3 &b, const NumericMatrix &R) {
	std::vector<double> pdiff(nloc);
	std::vector<double> psame(nloc);

	double puniq = A(i,ng); // mutation rate
	for (int l = 0; l < nloc; l++) {
		int human_allele = human[id][l];
		pdiff[l] = MAX(R(i,0) * b[i][l][human_allele],0.0);
		psame[l] = MAX(R(i,0) * b[i][l][human_allele] + R(i,1) * (1.0-A(i,ng)),0.0);
	}
	double lik = 0.0;
	for (int ii = 0; ii < ng; ii++) {								// Cycle through source of the clonal frame
		double mii = A(i,ii)/(1.0-A(i,ng));
		double l_ii = 0.0;
		for (int jj = 0; jj <nST[ii]; jj++) {
			double l_jj = mii;						//	Cycle through each ST from that source
			bool* HUMAN_UNIQUE = human_unique[id];
			bool* SAME = same[id][ii][jj];
			for(int l=0; l<nloc; l++, HUMAN_UNIQUE++, SAME++) {
				if(*HUMAN_UNIQUE) {						// new allele (allow some rounding error)
					l_jj *= puniq;
				}
				else if(*SAME) {						// previously observed and same as CF
					l_jj *= psame[l];
				}
				else {									// previously observed but different to CF
					l_jj *= pdiff[l];
				}
			}
			l_ii += l_jj * ABUN[ii][jj];
		}
		lik += l_ii / size[ii];
	}
	return lik;
}

void Cluster::precalc() {

	human_unique = Matrix<bool>(human.nrows(),nloc);
	beast_unique = Vector< Matrix<bool> >(ng);
	for (int i = 0; i < ng; i++)
	  beast_unique[i] = Matrix<bool>(nST[i],nloc);

	same = new bool***[human.nrows()];
	for(int i = 0; i < human.nrows(); i++) {
		same[i] = new bool**[ng];
		for(int ii = 0; ii < ng; ii++) {
			same[i][ii] = new bool*[nST[ii]];
			for(int jj = 0; jj < nST[ii]; jj++) {
				same[i][ii][jj] = new bool[nloc];
				for(int l = 0; l < nloc; l++) {
					same[i][ii][jj][l] = (human[i][l]==MLST[ii][jj][l]);
				}
			}
		}
		for(int l = 0; l < nloc; l++) {
			int human_allele = human[i][l];
		  human_unique[i][l] = (human_allele>=acount[ng][l].size()
                            || acount[ng][l][human_allele]==0);
		}
	}

	ksame = new bool****[ng];
	for(int i = 0; i < ng; i++) {
		ksame[i] = new bool***[nST[i]];
		for(int j = 0; j < nST[i]; j++) {
			ksame[i][j] = new bool**[ng];
			for(int ii = 0; ii < ng; ii++) {
				ksame[i][j][ii] = new bool*[nST[ii]];
				for(int jj = 0; jj < nST[ii]; jj++) {
					ksame[i][j][ii][jj] = new bool[nloc];
					for(int l = 0; l < nloc; l++) {
						ksame[i][j][ii][jj][l] = (MLST[i][j][l] == MLST[ii][jj][l]);
					}
				}
			}
			for(int l = 0; l < nloc; l++) {
				int allele = MLST[i][j][l];
				double num = acount[ng][l][allele] * (double)size[ng];
				beast_unique[i][j][l] = (num < 1.1);
			}
		}
	}
}

void append_traces(int iter, NumericMatrix &A, NumericMatrix &R, double lik, Matrix<double> &traces, int trace_row) {
  int col = 0;
  if (trace_row >= traces.nrows())
    return;
  traces[trace_row][col++] = iter;
  for (int i = 0; i < A.nrow(); i++) {
    for (int j = 0; j < A.ncol(); j++)
      traces[trace_row][col++] = A(i,j);
  }
  for (int i = 0; i < R.nrow(); i++)
    traces[trace_row][col++] = R(i,0);
  traces[trace_row][col++] = lik;
}

/* This version uses the clonal frame version of the likelihood */
void Cluster::mcmc6f(const double alpha, const double beta, const double gamma_, const int niter, const int thin, Random &ran) {
	precalc();

	int i,j;
	/* Initialize the Markov chain */
	int use = 0; int notuse = (int)!use;
	std::vector<double> BETA(ng+1,beta);			///<	Dirichlet hyperparameters of migration matrix (beta==1)
	NumericMatrix a(ng,ng+1);		///<	Reparametrised migration matrix. a[,ng]=mutation rate
	NumericMatrix A(ng,ng+1);		///<	Migration matrix M, M[ng] = mutation rates?
	const bool a_constraint = false;
	for(i=0;i<ng;i++) {
		while(true) {
			for(j=0;j<ng+1;j++) {
				a(i,j) = ran.gamma(1.,BETA[j]);
			}

			if(!a_constraint) break;
			double amax = a(i,0);
			for(j=1;j<ng;j++) if(a(i,j)>amax) amax = a(i,j);
			if(a(i,i)==amax) break;
		}
	}
	calc_A(a, A);	///< Does transformation to M
	Array3 b = calc_b(A);

	NumericMatrix r(ng,2);	///< Reparameterised per-group recombination rates
	NumericMatrix R(ng,2);  ///< R[grp,1:2] = r[grp,1:2]/sum(r[grp,1:2])
	for(i=0;i<ng;i++) {
		//r[i] = ran.beta(gamma_,gamma_);
		for(j=0;j<2;j++) {
			r(i,j) = ran.gamma(1.,gamma_);
		}
	}
	calc_R(r, R);

	/* Storage for likelihoods */
	double loglikelihood = known_source_loglik(A, b, R);

	/* Proposal probabilities */
	Vector<double> proprob(6,0.0);
	proprob[0] = 42./5.;							//	Update A: switching proposal
	proprob[1] = 42.;							//	Update A: log-normal proposal
	proprob[4] = 12./5.;							//	Update r: switching proposal
	proprob[5] = 12.;							//	Update r: log-normal proposal
	double tot = 0.0;
	for(i=0;i<proprob.size();i++) tot += proprob[i];
	for(i=0;i<proprob.size();i++) proprob[i] /= tot;

	double sigma_a = 0.5;							//	factor for normal proposal in MH change of a (case 1)
	double sigma_r = 0.5;							//	factor for normal proposal in MH change of r (case 5)

	/* Trace output matrix */
	evolution_traces.resize(niter/thin+1, 1+ng*(ng+1)+ng+1); // iter, A, r, loglik
	int trace_row = 0;
	append_traces(0, A, R, loglikelihood, evolution_traces, trace_row++);

	clock_t start = clock(), current;
	clock_t next = start + (clock_t)CLOCKS_PER_SEC;
	Rcpp::Rcout << "Done 0 of " << niter << " iterations";

	int iter;
	const int burnin = (int)floor((double)niter*.1);
	const int inc = MAX((int)floor((double)niter/100.),1);
	for(iter=0;iter<niter+burnin;iter++) {
		if(iter>=burnin && (iter-burnin)%inc==0) {

			/* Now dump the likelihoods for this iteration.
			   Weird that this has nothing to do with the thinning, which applies only
		     to the evolutionary parameters. */

			/* Compute likelihood of human isolate from each source */
		  Matrix<double> phi(human.nrows(), ng);
		  {
			  for(int h = 0; h < human.nrows(); h++) {
          // calculate the likelihood
          for (int j = 0; j < ng; j++) {
            phi[h][j] = log(likHi6(h, j, A, b, R));
          }
			  }
			}

			human_likelihoods[iter] = phi;
		}
		else {
			int move = multinom(proprob,ran);			//	random sweep for proposing moves
			switch(move) {
				case 0:	{// update A: switching proposal
					int popid = ran.discrete(0,ng-1);	// Choose the source for which to change the "mig" matrix
					int id1 = ran.discrete(0,ng);		// Principal elements of mig matrix to change
					int id2 = ran.discrete(0,ng-1);
					if(id2==id1) id2 = ng;
					if(a_constraint && (id1==popid || id2==popid)) break;
					// Need only update one row of a
					NumericMatrix::Row ar = a(popid,_);
					NumericVector ar_prop = ar;
					ar_prop[id1] = ar[id2];
					ar_prop[id2] = ar[id1];
					NumericMatrix A_prop(clone(A));
					A_prop(popid,_) = normalise(ar_prop);
					double logalpha = 0.0;
					// Prior ratio equals 1 because prior is symmetric
					// Hastings ratio equals 1 because proposal is symmetric
					// Likelihood ratio
					Array3 b_prop = calc_b(A_prop);
					double newloglik = known_source_loglik(A_prop, b_prop, R);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || ran.U() < exp(logalpha)) {	// accept
						ar = ar_prop;
					  A(popid,_) = A_prop(popid,_);
					  b  = b_prop;
						loglikelihood = newloglik;
					}
					else { // reject
					}
					break;
				}
				case 1:	{// update A: log-normal proposal
					int popid = ran.discrete(0,ng-1);	// Choose the source for which to change the "mig" matrix
					int id = ran.discrete(0,ng);		// Principal element of mig matrix to change
					// Need only update one row of a
					NumericMatrix::Row ar = a(popid,_);
					NumericVector ar_prop = ar;
					ar_prop[id] = exp(ran.normal(log(ar[id]), sigma_a));
					bool reject = false;
					if(a_constraint) {
						double ap_primemax = ar_prop[0];
						for(j=1;j<ng;j++) if(ar_prop[j]>ap_primemax) ap_primemax = ar_prop[j];
						if(ar_prop[popid]!=ap_primemax) reject = true;
					}
					if(reject) break;
					NumericMatrix A_prop(clone(A));
					A_prop(popid,_) = normalise(ar_prop);
					// Prior-Hastings ratio
					double logalpha = ar[id] - ar_prop[id];
					logalpha += beta * log(ar_prop[id]/ar[id]);
					// Likelihood ratio
					Array3 b_prop = calc_b(A_prop);
					double newloglik = known_source_loglik(A_prop, b_prop, R);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || ran.U() < exp(logalpha)) {	// accept
						ar = ar_prop;
					  A(popid,_) = A_prop(popid,_);
					  b  = b_prop;
						loglikelihood = newloglik;
					}
					else { // reject
					}
					break;
				}
				case 4: {// update r (switching move)
					int popid = ran.discrete(0,ng-1);
				  // Need only update one row of r
				  NumericMatrix::Row rr = r(popid,_);
				  NumericVector rr_prop = rr;
					rr_prop[0] = rr[1];
					rr_prop[1] = rr[0];
					NumericMatrix R_prop(clone(R));
					R_prop(popid,_) = normalise(rr_prop);
					double logalpha = 0.0;
					// Prior ratio equals 1 because prior is symmetric
					// Symmetric proposal so Hastings ratio equals 1
					// Likelihood ratio
					double newloglik = known_source_loglik(A, b, R_prop);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || ran.U() < exp(logalpha)) {	// accept
						rr = rr_prop;
						R(popid,_) = R_prop(popid,_);
						loglikelihood = newloglik;
					}
					else { // reject
					}
					break;
				}
				case 5:	{// update r (log-normal move)
					int popid = ran.discrete(0,ng-1);	// Choose the source for which to change the "rec" parameter
					int id = ran.discrete(0,1);			// Change one or other of the gamma components
					// Need only update one row of r
					NumericMatrix::Row rr = r(popid,_);
					NumericVector rr_prop = rr;
					rr_prop[id] = exp(ran.normal(log(rr[id]), sigma_r));
					NumericMatrix R_prop = clone(R);
					R_prop(popid,_) = normalise(rr_prop);
					// Prior-Hastings ratio
					double logalpha = rr[id] - rr_prop[id];
					logalpha += gamma_ * log(rr_prop[id]/rr[id]);
					// Likelihood ratio
					double newloglik = known_source_loglik(A, b, R_prop);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || ran.U() < exp(logalpha)) {	// accept
						rr = rr_prop;
						R(popid, _) = R_prop(popid, _);
						loglikelihood = newloglik;
					}
					else { // reject
					}
					break;
				}
				default: {
					error("Move not recognised");
				}
			}
		}

		/* output thinned traces of island model fit */
		if(iter >= burnin && (iter+1)%thin==0)
		  append_traces(iter+1, A, R, loglikelihood, evolution_traces, trace_row++);

		if((current=clock())>next) {
		  Rcpp::Rcout << "\rDone " << (iter+1) << " of " << niter+burnin << " iterations in " << (double)(current-start)/CLOCKS_PER_SEC << " s " << std::flush;
			next = current + (clock_t)CLOCKS_PER_SEC;
		}
	}
	Rcpp::Rcout << std::endl;
}

// helper stuff below here

// Assumes A is correctly sized
void Cluster::calc_A(NumericMatrix &a, NumericMatrix &A) {
  for (int i = 0; i < a.nrow(); i++) {
    A(i,_) = normalise(a(i,_));
  }
}

// Assumes R is correctly sized
void Cluster::calc_R(NumericMatrix &r, NumericMatrix &R) {
  for (int i = 0; i < r.nrow(); i++) {
    R(i,_) = normalise(r(i,_));
  }
}

NumericVector Cluster::normalise(const NumericVector &x) {
  return x / sum(x);
}

Cluster::Array3 Cluster::calc_b(const NumericMatrix &A) {
  Array3 b(ng);
  for(int i = 0; i < ng; i++) {
    b[i].resize(nloc);
    for(int j = 0; j < nloc; j++) {
      b[i][j].resize(acount[i][j].size());
      for(int k = 0; k < acount[i][j].size(); k++) {
        b[i][j][k] = 0.0;
        for(int l = 0; l < ng; l++) {
          b[i][j][k] += acount[l][j][k] * A(i,l);
          //					bk[i][j][k] += (acount[l][j][k]*(double)size[l]-1.0)/(double)(size[l]-1) * a[i][l];
        }
      }
    }
  }
  return b;
}

int Cluster::multinom(Vector<double> &p, Random &ran) {
  double U = ran.U();
  int i;
  for(i=0;i<p.size();i++) {
    if((U-=p[i])<=0.0) break;
  }
  if(U>0.0) error("Problem in multinom");
  return i;
}

void Cluster::initialise(Matrix<int> &isolate) {
  init = true;

  Rcpp::Rcout << "Running initialise..." << std::endl;

  /* Format assumed is ST <genes> SOURCE */
  nloc = isolate.ncols() - 2;
  ng   = 0;
  for (int i = 0; i < isolate.nrows(); i++) {
    const int sc = isolate[i][nloc+1];
    if (sc > ng)
      ng = sc;
  }

  // find maximum ST, maximum allele for each locus, the total in each source group, humans, and across all sources
  int maxST = -1;			///< maximum ST
  Vector<int> maxallele(nloc,-1);	///< maximum allele number for each locus
  size = Vector<int>(ng+1,0);	///< size[0:(ng-1)] is total in each source group, size[ng] is total in source groups
  int nhuman = 0; 		///< number of human isolates
  for (int i = 0; i < isolate.nrows(); i++) {
    if (isolate[i][0] > maxST)
      maxST = isolate[i][0];
    for (int j = 1; j <= nloc; j++) {
      if (isolate[i][j] > maxallele[j-1])
        maxallele[j-1] = isolate[i][j];
    }
    if (isolate[i][nloc+1] > 0) {
      ++size[isolate[i][nloc+1]-1];
      ++size[ng];
    }
    else
      ++nhuman;
  }
  Rcpp::Rcout << "Maximum ST is " << maxST << " total isolates is " << size[ng] << std::endl;

  // map STs to their numbers
  Matrix<int> aprofile(MIN(maxST,isolate.nrows()),nloc+1,-1);
  Vector<int> STwhere(maxST+1,-1);
  int NST = 0;
  for (int i=0; i < isolate.nrows(); i++) {
    const int lab = isolate[i][0];
    // have we seen this one before?
    if (STwhere[lab] == -1) {
      // nope - add it
      for (int j = 0; j < nloc; j++)
        aprofile[NST][j] = isolate[i][j+1];
      STwhere[lab] = NST;
      ++NST;
    }
  }
  Rcpp::Rcout << "Number of STs is " << NST << " number of humans is " << nhuman << " Number of loci is " << nloc << " Number of groups is " << ng << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Create the matrix of human isolates
  human.resize(nhuman,nloc);
  int ih = 0;
  for (int i=0; i < isolate.nrows(); i++) {
    if (isolate[i][nloc+1] == 0) {
      for (int j = 0; j < nloc; j++) {
        human[ih][j] = isolate[i][j+1];
      }
      ++ih;
    }
  }

  // Right, down here we should have everything we need:
  //  humans  - Matrix of MLST genes, complete with duplicates
  //  STwhere - map from ST to number (0..NST-1)
  //  isolate - actual isolate data of the form [ST MLST_genes Group]
  //  NST     - number of STs

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Count the number of non-human isolates in each group
  /* Calculate the size and number of unique STs in each group */
  Matrix<int> niso(NST,ng+1,0);	///< number of each ST in each group, niso[,ng] is number of each ST in all groups
  size = Vector<int>(ng+1,0);	///< number of STs in each group, size[ng] is number of STs
  nST = Vector<int>(ng+1,0);	///< number of unique STs in each group, nST[ng] is number of unique STs
  for (int i = 0; i < isolate.nrows(); i++) {
    const int lab = isolate[i][0];
    const int gp  = isolate[i][nloc+1]-1;
    const int wh  = STwhere[lab];
    if (gp >= 0) {
      if (niso[wh][gp] == 0) ++nST[gp];
      if (niso[wh][ng] == 0) ++nST[ng];
      ++niso[wh][gp];
      ++niso[wh][ng];
      ++size[gp];
      ++size[ng];
    }
  }

  /* Allocate memory for allele counts */
  acount.resize(ng+1,nloc);	///< counts of alleles in each group at each loci
  for (int i = 0; i < ng+1; i++) {
    for (int j = 0; j < nloc; j++) {
      acount[i][j].resize(maxallele[j]+1);
      for (int k = 0; k <= maxallele[j]; k++)
        acount[i][j][k] = 0;
    }
  }
  /* Record the allelic profile for each unique ST in each group */
  MLST.resize(ng);		///< MLST[group,unique_st,loci] -> MLST[group] is profile of each unique ST
  FREQ.resize(ng);		///< FREQ[group,unique_st] proportion of unique_st in group
  ABUN.resize(ng);		///< ABUN[group,unique_st] number of unique_st in group
  for (int i = 0; i < ng; i++) {
    MLST[i].resize(nST[i],nloc);
    FREQ[i].resize(nST[i]);
    ABUN[i].resize(nST[i]);
  }
  Vector<int> ix(ng,0);		///< counter for each group - the allelic profile is copied in separately for each group
  for (int i = 0; i < NST; i++) { // for each unique ST
    for (int sc = 0; sc < ng; sc++) { // for each group
      const int ct = niso[i][sc]; // how many of this ST is in this group?
      if (ct>0) {
        for (int j = 0; j < nloc; j++) { // copy across the MLST profile
          MLST[sc][ix[sc]][j] = aprofile[i][j];
        }
        FREQ[sc][ix[sc]] = (double)ct/(double)size[sc];
        ABUN[sc][ix[sc]] = (double)ct;
        for(int j = 0; j < nloc; j++) { // for each loci
          const int allele = aprofile[i][j];			// allele
          acount[sc][j][allele] += (double)ct/(double)size[sc];	// allele frequency
          //acount[ng][j][allele] += ct * (weight = size[ng]/ng/size[sc]) /size[ng]
          /* weighted */// acount[ng][j][allele] += (double)ct/(double)ng/(double)size[sc];
          /* unweighted */ acount[ng][j][allele] += (double)ct/(double)size[ng];
        }
        ++ix[sc];
      }
    }
  }

  nalleles = Matrix<int>(ng+1,nloc,0);	///< number of alleles for each group and each loci. nalleles[ng] is total
  for (int i = 0; i <= ng; i++) {
    for (int j = 0; j < nloc; j++) {
      for (int k = 0; k <= maxallele[j]; k++) {
        if (acount[i][j][k]>0) {
          ++nalleles[i][j];
        }
      }
    }
  }
}

