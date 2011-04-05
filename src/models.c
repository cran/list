#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

void ictregBinom(int *Y,             /* outcome vector */
		 int *J,             /* # of control items */
		 int *n_samp,        /* sample size */
		 int *n_draws,       /* # of MCMC draws */
		 int *treat,         /* treatment indicator vector: 0 or 1 */
		 double *Xall,       /* covariates in a vector form */
		 double *delta,      /* coefs for sensitive item */
		 double *psi,        /* coefs for control items */ 
		 int *n_cov,         /* # of covariates */
		 double *delta0,     /* prior mean for delta */
		 double *psi0,       /* prior mean for psi */
		 double *A0deltaAll, /* prior precision for delta */
		 double *A0psiAll,   /* prior precision for psi */
		 double *deltaVar,   /* proposal variance for delta */
		 double *psiVar,     /* proposal variance for psi */
		 int *unconst,       /* is this unconstrained model? */
		 int *burnin,        /* number of burnins */
		 int *keep,          /* keep every *th draw */
		 int *verbose,       /* want to print progress? */
		 double *allresults  /* storage for results */
		 ) {

  int i, j, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int itempZstar0, itempZstar1;
  int *deltaCounter = intArray(1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);    /* acceptance ratio for psi */
  int *psi1Counter = intArray(1);   /* acceptance ratio for psi1 */
  double dtemp1, dtemp2;

  /** get random seed **/
  GetRNGstate();
  /* additional parameters */
  double *psi1 = doubleArray(*n_cov);
  double *psi01 = doubleArray(*n_cov);
  
  if (*unconst) {
    for (j = 0; j < *n_cov; j++) {
      psi1[j] = psi[*n_cov + j];
      psi01[j] = psi0[*n_cov + j];;
    }
  }

  /** Data **/
  int *Zstar = intArray(*n_samp);
  int *Y0 = intArray(*n_samp);
  int *Y0temp0 = intArray(*n_samp);
  int *Y0temp1 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, *n_cov); 
  double **Xtemp0 = doubleMatrix(*n_samp, *n_cov); 
  double **Xtemp1 = doubleMatrix(*n_samp, *n_cov); 
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double *Xpsi1 = doubleArray(*n_samp);
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  /** Prior **/
  double **A0delta = doubleMatrix(*n_cov, *n_cov);
  double **A0psi = doubleMatrix(*n_cov, *n_cov);
  double **A0psi1 = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0delta[i][j] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  if (*unconst) {
    itemp = 0;
    for (j = 0; j < *n_cov; j++)
      for (i = 0; i < *n_cov; i++)
	A0psi1[i][j] = A0psiAll[itemp++];
  }

  /** Proposal precisoin **/
  double **deltaPro = doubleMatrix(*n_cov, *n_cov);
  double **psiPro = doubleMatrix(*n_cov, *n_cov);
  double **psi1Pro = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      deltaPro[i][j] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      psiPro[i][j] = psiVar[itemp++];

  if (*unconst) {
    itemp = 0;
    for (j = 0; j < *n_cov; j++)
      for (i = 0; i < *n_cov; i++)
	psi1Pro[i][j] = psiVar[itemp++];
  }

  /*
  PdoubleMatrix(A0delta, *n_cov, *n_cov);
  PdoubleMatrix(A0psi, *n_cov, *n_cov);
  if (*unconst) {
    PdoubleMatrix(A0psi1, *n_cov, *n_cov);
  }
  PdoubleMatrix(deltaPro, *n_cov, *n_cov);
  PdoubleMatrix(psiPro, *n_cov, *n_cov); 
  if (*unconst) {
    PdoubleMatrix(psi1Pro, *n_cov, *n_cov);
  } 
  */

  /** MCMC **/
  itemp = 0; deltaCounter[0] = 0; psiCounter[0] = 0; psi1Counter[0] = 0;
  if (*verbose) {
    Rprintf("Starting posterior sampling...\n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    itempZstar0 = 0; itempZstar1 = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar */
      if ((treat[i] == 1) && (Y[i] == (*J+1))) {
	Zstar[i] = 1;
	Y0[i] = *J;
      } else if ((treat[i] == 1) && (Y[i] == 0)) {
	Zstar[i] = 0;
	Y0[i] = 0;
      } else {
	Xdelta[i] = 0;  Xpsi[i] = 0;  Xpsi1[i] = 0;
	for (j = 0; j < *n_cov; j++) {
	  Xdelta[i] += X[i][j]*delta[j];
	  Xpsi[i] += X[i][j]*psi[j];
	  if (*unconst) {
	    Xpsi1[i] += X[i][j]*psi1[j];
	  }
	}
	if (*unconst) {
	  dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
		       dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi1[i])), 1));
	} else {
	  dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
		       dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	}
	dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	  Zstar[i] = 1;
	  Y0[i] = Y[i] - treat[i];
	} else { 
	  Zstar[i] = 0;
	  Y0[i] = Y[i];
	}
	if (*unconst) {
	  if (Zstar[i] == 1) {
	    Y0temp0[itempZstar0] = Y0[i];
	    for (j = 0; j < *n_cov; j++) {
	      Xtemp0[itempZstar0][j] = X[i][j];
	    }
	    itempZstar0++;
	  } else {
	    Y0temp1[itempZstar1] = Y0[i];
	    for (j = 0; j < *n_cov; j++) {
	      Xtemp1[itempZstar1][j] = X[i][j];
	    }
	    itempZstar1++;
	  }	    
	}
      }
    }
    
    /* Sample delta */
    BinomLogitMetro(Zstar, X, delta, *n_samp, 1, *n_cov, delta0, A0delta, deltaPro, 1, deltaCounter);
    
    /* Sample psi */
    if (*unconst) {
      BinomLogitMetro(Y0temp0, Xtemp0, psi, itempZstar0, *J, *n_cov, psi0, A0psi, psiPro, 1, psiCounter);
      BinomLogitMetro(Y0temp1, Xtemp1, psi1, itempZstar1, *J, *n_cov, psi01, A0psi1, psi1Pro, 1, psi1Counter);
    } else {
      BinomLogitMetro(Y0, X, psi, *n_samp, *J, *n_cov, psi0, A0psi, psiPro, 1, psiCounter);
    }

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = delta[j];
	}      
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = psi[j];
	}
	if (*unconst) {
	  for (j = 0; j < *n_cov; j++) {
	    allresults[itemp++] = psi1[j];
	  }
	}
	allresults[itemp++] = ((double) *deltaCounter / (double) (main_loop + 1));
	allresults[itemp++] = ((double) *psiCounter / (double) (main_loop + 1));
	if (*unconst) {
	  allresults[itemp++] = ((double) *psi1Counter / (double) (main_loop + 1));
	}
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(Zstar);
  free(Y0);
  free(Y0temp0);
  free(Y0temp1);
  free(Xdelta);
  free(Xpsi);
  free(Xpsi1);
  free(psi1);
  free(psi01);
  free(deltaCounter);
  free(psiCounter);
  free(psi1Counter);
  FreeMatrix(X, *n_samp);
  FreeMatrix(Xtemp0, *n_samp);
  FreeMatrix(Xtemp1, *n_samp);
  FreeMatrix(A0delta, *n_cov);
  FreeMatrix(A0psi, *n_cov);
  FreeMatrix(A0psi1, *n_cov);
  FreeMatrix(deltaPro, *n_cov);
  FreeMatrix(psiPro, *n_cov);
  FreeMatrix(psi1Pro, *n_cov);
}


/*** 
   A Random Walk Metroplis Sampler for Binomial Logistic Regression 
   with Independent Normal Prior
   
   proposal distribution is the univariate normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.
***/

void BinomLogitMetro(int *Y,        /* outcome variable: 0, 1, ..., J */
		     double **X,    /* (N x K) covariate matrix */
		     double *beta,  /* K coefficients */
		     int n_samp,    /* # of obs */
		     int n_size,    /* # of size, J */
		     int n_cov,     /* # of covariates, K */
		     double *beta0, /* K prior mean vector */
		     double **A0,   /* (K x K) prior precision */
		     double **Var, /* K proposal precision */
		     int n_gen,     /* # of MCMC draws */
		     int *counter   /* # of acceptance */
		     ) {
  
  int i, j, main_loop;
  double numer, denom, Xbeta, Xprop;
  double *prop = doubleArray(n_cov);

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
    /** Sample from the proposal distribution **/
    rMVN(prop, beta, Var, n_cov);
    
    /** Calculating the ratio (log scale) **/
    /* prior */
    numer = dMVN(prop, beta0, A0, n_cov, 1);
    denom = dMVN(beta, beta0, A0, n_cov, 1);   
    
    /* likelihood */
    for (i = 0; i < n_samp; i++) {
      Xbeta = 0;
      Xprop = 0;
      for (j = 0; j < n_cov; j++) {
	Xbeta += X[i][j]*beta[j];
	Xprop += X[i][j]*prop[j];
      }
      denom += dbinom(Y[i], n_size, 1 / (1 + exp(-Xbeta)), 1); 
      numer += dbinom(Y[i], n_size, 1 / (1 + exp(-Xprop)), 1); 
    }
      
    /** Rejection **/
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      counter[0]++;
      for (j = 0; j < n_cov; j++) {
	beta[j] = prop[j];
      }
    }
  }
  
  free(prop);
} /* end of BinomlogitMetro */


