#include "asymmetric_island.h"

#include <sstream> // for the stream stuff
#include <iostream>

#include <Rcpp.h>

#include "logdouble.h"

using namespace Rcpp;

double Island::known_source_loglik(const NumericMatrix &A, const NumericArray3 &b, const NumericMatrix &R) {
	logdouble lik = 1.0;
	/* Cycle through each unique ST in each group, taking account of abundance of the STs */
	for (int i = 0; i < ng; i++) {
		double punique = A(i,ng);
		for (int j = 0; j < nST[i]; j++) {
			double ncopiesj = ABUN[i][j];
			std::vector<double> psame(nloc);
			std::vector<double> pdiff(nloc);
			for (int l = 0; l < nloc; l++) {
				int allele = MLST[i](j, l);
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
			std::vector<logdouble> l_j(ng);
			for (int ii = 0; ii < ng; ii++) {						//	Cycle through source of the clonal frame
				double mii = A(i,ii)/(1.0-A(i,ng));
				std::vector<logdouble> l_ii(nST[ii]);      // allocate the vector
				for (int jj = 0; jj < nST[ii]; jj++) {				//	Cycle through each ST from that source
					double ncopiesjj = (i==ii && j==jj) ? ABUN[ii][jj]-std::min(ABUN[ii][jj],1.0)
						: ABUN[ii][jj]; // NOTE: This can be zero.
					logdouble l_jj = mii;
					LogicalMatrix::Row BEAST_UNIQUE = beast_unique[i](j, _);
					bool *SAME = ksame[i][j][ii][jj];
					for (int l = 0; l < nloc; l++, SAME++) {
						if (BEAST_UNIQUE[l]) {				// new allele (allow some rounding error)
							l_jj *= punique;
						}
						else if(*SAME) {				// previously observed and same as CF
							l_jj *= psame[l];
						}
						else {							// previously observed but different to CF
							l_jj *= pdiff[l];
						}
					}
					l_ii[jj] = l_jj * ncopiesjj;
				}
				l_j[ii] = sum(l_ii) / size[ii];
			}
			lik *= sum(l_j)^ncopiesj;
		}
	}
	return lik.log();
}

double Island::likHi6(const int id, const int i, const NumericMatrix &A, const NumericArray3 &b, const NumericMatrix &R) {
	std::vector<double> pdiff(nloc);
	std::vector<double> psame(nloc);

	double puniq = A(i,ng); // mutation rate
	for (int l = 0; l < nloc; l++) {
		int human_allele = human(id, l);
		pdiff[l] = std::max(R(i,0) * b[i][l][human_allele],0.0);
		psame[l] = std::max(R(i,0) * b[i][l][human_allele] + R(i,1) * (1.0-A(i,ng)),0.0);
	}
	double lik = 0.0;
	for (int ii = 0; ii < ng; ii++) {								// Cycle through source of the clonal frame
		double mii = A(i,ii)/(1.0-A(i,ng));
		double l_ii = 0.0;
		for (int jj = 0; jj <nST[ii]; jj++) {
			double l_jj = mii;						//	Cycle through each ST from that source
			LogicalMatrix::Row HUMAN_UNIQUE = human_unique(id, _);
			bool* SAME = same[id][ii][jj];
			for(int l=0; l<nloc; l++, SAME++) {
				if (HUMAN_UNIQUE[l]) {						// new allele (allow some rounding error)
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

void Island::precalc() {
  // TODO: Much of this can probably be vectorised

  // NOTE: This just plain doesn't work for those human STs that are _identical_ to
  //       beast STs. i.e. we can't attribute beast STs using this scenario, as
  //       if we're pretending a beast ST is human, then human_unique will always
  //       be false (as we would have seen all alleles before). Similarly, same
  //       will always be true, as we would have seen it before.
  //       So the ST distribution we derive only works for those types we haven't
  //       seen before, or are genuinely identical between beast and human.
  //
  //       The correct thing to do if we wish to also get the ST distribution for
  //       non-human STs is compute it using beast_unique etc.

	human_unique = LogicalMatrix(human.nrow(), nloc);
  for(int i = 0; i < human.nrow(); i++) {
    for(int l = 0; l < nloc; l++) {
      int human_allele = human(i, l);
      human_unique(i,l) = (human_allele>=(int)acount[ng][l].size()
                             || acount[ng][l][human_allele]==0);
    }
  }

  beast_unique.resize(ng);
	for (int i = 0; i < ng; i++) {
	  beast_unique[i] = LogicalMatrix(nST[i], nloc);
	  for(int j = 0; j < nST[i]; j++) {
	    for(int l = 0; l < nloc; l++) {
	      int allele = MLST[i](j,l);
	      double num = acount[ng][l][allele] * (double)size[ng];
	      beast_unique[i](j,l) = (num < 1.1);
	    }
	  }
	}

	same = new bool***[human.nrow()];
	for(int i = 0; i < human.nrow(); i++) {
		same[i] = new bool**[ng];
		for(int ii = 0; ii < ng; ii++) {
			same[i][ii] = new bool*[nST[ii]];
			for(int jj = 0; jj < nST[ii]; jj++) {
				same[i][ii][jj] = new bool[nloc];
				for(int l = 0; l < nloc; l++) {
					same[i][ii][jj][l] = (human(i, l) == MLST[ii](jj,l));
				}
			}
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
						ksame[i][j][ii][jj][l] = (MLST[i](j,l) == MLST[ii](jj,l));
					}
				}
			}
		}
	}
}

void append_traces(int iter, NumericMatrix &A, NumericMatrix &R, double lik, NumericMatrix &traces, int trace_row) {
  int col = 0;
  if (trace_row >= traces.nrow())
    return;
  NumericMatrix::Row row = traces(trace_row,_);
  row[col++] = iter;
  for (int i = 0; i < A.nrow(); i++) {
    for (int j = 0; j < A.ncol(); j++)
      row[col++] = A(i,j);
  }
  for (int i = 0; i < R.nrow(); i++)
    row[col++] = R(i,0);
  row[col++] = lik;
}

/* This version uses the clonal frame version of the likelihood */
void Island::mcmc6f(const double beta, const double gamma_, const int samples, const int burnin, const int thin) {
	precalc();

  /* Initialize the random number generator */
  RNGScope scope;

	/* Initialize the Markov chain */
	NumericMatrix a(ng,ng+1);		///<	Reparametrised migration matrix. a[,ng]=mutation rate
	NumericMatrix A(ng,ng+1);		///<	Migration matrix M, M[ng] = mutation rates?
	for (int i = 0; i < ng; i++) {
	  a(i,_) = rgamma(ng+1, beta, 1.0);
	}
	calc_A(a, A);	///< Does transformation to M
	NumericArray3 b = calc_b(A);

	NumericMatrix r(ng,2);	///< Reparameterised per-group recombination rates
	NumericMatrix R(ng,2);  ///< R[grp,1:2] = r[grp,1:2]/sum(r[grp,1:2])
	for (int i = 0; i < ng; i++) {
		r(i,_) = rgamma(2, gamma_, 1.0);
	}
	calc_R(r, R);

	/* Storage for likelihoods */
	double loglikelihood = known_source_loglik(A, b, R);

	/* Proposal probabilities */
	NumericVector proprob = NumericVector::create(42./5.,   //	Update A: switching proposal
                                                42.,      //	Update A: log-normal proposal
                                                0.0, 0.0,
                                                12./5.,   //	Update r: switching proposal
                                                12.);     //	Update r: log-normal proposal
  proprob = proprob / sum(proprob);

	double sigma_a = 0.5;							//	factor for normal proposal in MH change of a (case 1)
	double sigma_r = 0.5;							//	factor for normal proposal in MH change of r (case 5)

	/* Trace output matrix */
	evolution_traces = NumericMatrix(samples+1, 1+ng*(ng+1)+ng+1); // iter, A, r, loglik
	int trace_row = 0;
	append_traces(0, A, R, loglikelihood, evolution_traces, trace_row++);

//	clock_t start = clock(), current;
//	clock_t next = start + (clock_t)CLOCKS_PER_SEC;
//	Rcout << "Done 0 of " << niter << " iterations";

	for (int iter = 0; iter < (samples+burnin)*thin; iter++) {
		if ((iter+1) % thin == 0 && iter/thin >= burnin) {

		  /* JM Changed: Thinning is now the same - we just want given number of samples */
			/* Dump the likelihoods for this iteration.
			   Weird that this has nothing to do with the thinning, which applies only
		     to the evolutionary parameters. */

			/* Compute likelihood of human isolate from each source */
		  NumericMatrix phi(human.nrow(), ng);
		  for(int h = 0; h < human.nrow(); h++) {
        // calculate the likelihood
        for (int j = 0; j < ng; j++) {
          phi(h,j) = likHi6(h, j, A, b, R);
        }
		  }

		  std::stringstream s; s << iter;
			human_likelihoods[s.str()] = phi;
		}
		else {
			int move = multinom(proprob);			//	random sweep for proposing moves
			switch(move) {
				case 0:	{// update A: switching proposal
					int popid = sample(ng);	// Choose the source for which to change the "mig" matrix
					int id1 = sample(ng+1);		// Principal elements of mig matrix to change
					int id2 = sample(ng);
					if(id2==id1) id2 = ng;
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
					NumericArray3 b_prop = calc_b(A_prop);
					double newloglik = known_source_loglik(A_prop, b_prop, R);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || R::runif(0, 1) < exp(logalpha)) {	// accept
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
					int popid = sample(ng);	// Choose the source for which to change the "mig" matrix
					int id = sample(ng+1);		// Principal element of mig matrix to change
					// Need only update one row of a
					NumericMatrix::Row ar = a(popid,_);
					NumericVector ar_prop = ar;
					ar_prop[id] = R::rlnorm(log(ar[id]), sigma_a);

					NumericMatrix A_prop(clone(A));
					A_prop(popid,_) = normalise(ar_prop);
					// Prior-Hastings ratio
					double logalpha = ar[id] - ar_prop[id];
					logalpha += beta * log(ar_prop[id]/ar[id]);
					// Likelihood ratio
					NumericArray3 b_prop = calc_b(A_prop);
					double newloglik = known_source_loglik(A_prop, b_prop, R);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || R::runif(0, 1) < exp(logalpha)) {	// accept
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
					int popid = sample(ng);
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
					if (logalpha >= 0.0 || R::runif(0, 1) < exp(logalpha)) {	// accept
						rr = rr_prop;
						R(popid,_) = R_prop(popid,_);
						loglikelihood = newloglik;
					}
					else { // reject
					}
					break;
				}
				case 5:	{// update r (log-normal move)
					int popid = sample(ng);	// Choose the source for which to change the "rec" parameter
					int id = sample(2);			// Change one or other of the gamma components
					// Need only update one row of r
					NumericMatrix::Row rr = r(popid,_);
					NumericVector rr_prop = rr;
					rr_prop[id] = R::rlnorm(log(rr[id]), sigma_r);
					NumericMatrix R_prop = clone(R);
					R_prop(popid,_) = normalise(rr_prop);
					// Prior-Hastings ratio
					double logalpha = rr[id] - rr_prop[id];
					logalpha += gamma_ * log(rr_prop[id]/rr[id]);
					// Likelihood ratio
					double newloglik = known_source_loglik(A, b, R_prop);

					logalpha += newloglik - loglikelihood;
					if (logalpha >= 0.0 || R::runif(0, 1) < exp(logalpha)) {	// accept
						rr = rr_prop;
						R(popid, _) = R_prop(popid, _);
						loglikelihood = newloglik;
					}
					else { // reject
					}
					break;
				}
				default: {
				  stop("Move not recognised");
				}
			}
		}

		/* output thinned traces of island model fit */
		if((iter+1)%thin==0 && iter/thin >= burnin)
		  append_traces(iter+1, A, R, loglikelihood, evolution_traces, trace_row++);

//		if((current=clock())>next) {
//		  Rcout << "\rDone " << (iter+1) << " of " << niter+burnin << " iterations in " << (double)(current-start)/CLOCKS_PER_SEC << " s " << std::flush;
//			next = current + (clock_t)CLOCKS_PER_SEC;
//		}
	}
	Rcout << std::endl;
}

// helper stuff below here

// Assumes A is correctly sized
void Island::calc_A(NumericMatrix &a, NumericMatrix &A) {
  for (int i = 0; i < a.nrow(); i++) {
    A(i,_) = normalise(a(i,_));
  }
}

// Assumes R is correctly sized
void Island::calc_R(NumericMatrix &r, NumericMatrix &R) {
  for (int i = 0; i < r.nrow(); i++) {
    R(i,_) = normalise(r(i,_));
  }
}

NumericVector Island::normalise(const NumericVector &x) {
  return x / sum(x);
}

Island::NumericArray3 Island::calc_b(const NumericMatrix &A) {
  NumericArray3 b(ng);
  for(int i = 0; i < ng; i++) {
    b[i].resize(nloc);
    for(int j = 0; j < nloc; j++) {
      b[i][j].resize(acount[i][j].size());
      for(size_t k = 0; k < acount[i][j].size(); k++) {
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

int Island::sample(int n) {
  double U = R::runif(0, 1);
  return U * n;
}

int Island::multinom(const NumericVector &p) {
  double U = R::runif(0, 1);
  for (int i = 0; i < p.size(); i++) {
    if ((U-=p[i]) <= 0.0)
      return i;
  }
  stop("Problem in multinom");
  return 0;
}

void Island::initialise(IntegerMatrix isolate) {
  init = true;

  /* Format assumed is ST <genes> SOURCE */
  nloc = isolate.ncol() - 2;
  IntegerMatrix::Column st     = isolate(_, 0);
  IntegerMatrix loci           = isolate(_, Range(1, nloc));
  IntegerMatrix::Column source = isolate(_, nloc+1);
  ng   = max(source);

  // find maximum ST, maximum allele for each locus, the total in each source group, humans, and across all sources
  int maxST = max(st);			      ///< maximum ST
  IntegerVector maxallele(nloc);	///< maximum allele number for each locus
  for (int j = 0; j < nloc; j++) {
    maxallele[j] = max(loci(_, j));
  }
  // count the number of isolates
  size = IntegerVector(ng);	///< size[0:(ng-1)] is total in each source group
  for (int i = 0; i < ng; i++) {
    size[i] = std::count(source.begin(), source.end(), i+1);
  }
  size.push_back(sum(size)); ///< size[ng] is total in source groups
  int nhuman = std::count(source.begin(), source.end(), 0);

  Rcout << "Maximum ST is " << maxST << " total isolates is " << size[ng] << std::endl;

  // map STs to their numbers
  std::vector<IntegerVector> aprofile;
  IntegerVector STwhere(maxST+1, -1); ///< Map from ST to aprofile row number.
  for (int i=0; i < isolate.nrow(); i++) {
    const int lab = st[i];
    // have we seen this one before?
    if (STwhere[lab] == -1) {
      // nope - add it
      STwhere[lab] = aprofile.size();
      aprofile.push_back(loci(i, _));
    }
  }
  int NST = aprofile.size();
  Rcout << "Number of STs is " << NST << " number of humans is " << nhuman << " Number of loci is " << nloc << " Number of groups is " << ng << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Create the matrix of human isolates
  human = IntegerMatrix(nhuman, nloc);
  int ih = 0;
  for (int i=0; i < loci.nrow(); i++) {
    if (source[i] == 0) {
      human(ih++, _) = loci(i, _);
    }
  }

  // Right, down here we should have everything we need:
  //  humans  - Matrix of MLST genes, complete with duplicates
  //  STwhere - map from ST to number (0..NST-1)
  //  isolate - actual isolate data of the form [ST MLST_genes Group]
  //  NST     - number of STs

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Count the number of non-human isolates in each group
  /* Calculate the number of unique STs in each group */
  IntegerMatrix niso(NST, ng+1);	///< number of each ST in each group, niso[,ng] is number of each ST in all groups
  nST = IntegerVector(ng+1);	    ///< number of unique STs in each group, nST[ng] is number of unique STs
  for (int i = 0; i < isolate.nrow(); i++) {
    const int gp  = source[i] - 1;
    const int wh  = STwhere[st[i]];
    if (gp >= 0) {
      if (niso(wh,gp) == 0) ++nST[gp];
      if (niso(wh,ng) == 0) ++nST[ng];
      ++niso(wh,gp);
      ++niso(wh,ng);
    }
  }

  /* Allocate memory for allele counts */
  acount.resize(ng+1);	///< counts of alleles in each group at each loci
  for (int i = 0; i < ng+1; i++) {
    acount[i].resize(nloc);
    for (int j = 0; j < nloc; j++) {
      acount[i][j].resize(maxallele[j]+1);
    }
  }
  /* Record the allelic profile for each unique ST in each group */
  MLST.resize(ng);		///< MLST[group,unique_st,loci] -> MLST[group] is profile of each unique ST
  FREQ.resize(ng);		///< FREQ[group,unique_st] proportion of unique_st in group
  ABUN.resize(ng);		///< ABUN[group,unique_st] number of unique_st in group
  for (int i = 0; i < ng; i++) {
    MLST[i] = IntegerMatrix(nST[i],nloc);
    FREQ[i].resize(nST[i]);
    ABUN[i].resize(nST[i]);
  }
  IntegerVector ix(ng,0);		///< counter for each group - the allelic profile is copied in separately for each group
  for (int i = 0; i < NST; i++) { // for each unique ST
    for (int sc = 0; sc < ng; sc++) { // for each group
      const int ct = niso(i, sc); // how many of this ST is in this group?
      if (ct > 0) {
        MLST[sc](ix[sc],_) = aprofile[i];
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
}

