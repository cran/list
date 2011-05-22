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

/** 
  Item Count Technique Binomial Regression for the Standard Design
**/

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

  int i, j, n_covp, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int itempZstar0, itempZstar1;
  int *deltaCounter = intArray(1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);    /* acceptance ratio for psi */
  int *psi1Counter = intArray(1);   /* acceptance ratio for psi1 */
  double dtemp1, dtemp2;

  /* intercept only model */
  if (*unconst == 2) {
    n_covp = *n_cov + 1;
  } else {
    n_covp = *n_cov;
  }

  /** get random seed **/
  GetRNGstate();
  /* additional parameters */
  double *psi1 = doubleArray(n_covp);
  double *psi01 = doubleArray(n_covp);
  
  if (*unconst == 1) {
    for (j = 0; j < n_covp; j++) {
      psi1[j] = psi[n_covp + j];
      psi01[j] = psi0[n_covp + j];;
    }
  }

  /** Data **/
  int *Zstar = intArray(*n_samp);
  int *Y0 = intArray(*n_samp);
  int *Y0temp0 = intArray(*n_samp);
  int *Y0temp1 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_covp); 
  double **Xtemp0 = doubleMatrix(*n_samp, n_covp); 
  double **Xtemp1 = doubleMatrix(*n_samp, n_covp); 
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double *Xpsi1 = doubleArray(*n_samp);
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  /** Prior **/
  double **A0delta = doubleMatrix(*n_cov, *n_cov);
  double **A0psi = doubleMatrix(n_covp, n_covp);
  double **A0psi1 = doubleMatrix(n_covp, n_covp);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0delta[i][j] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < n_covp; j++)
    for (i = 0; i < n_covp; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  if (*unconst == 1) {
    itemp = 0;
    for (j = 0; j < n_covp; j++)
      for (i = 0; i < n_covp; i++)
	A0psi1[i][j] = A0psiAll[itemp++];
  }

  /** Proposal precision **/
  double **deltaPro = doubleMatrix(*n_cov, *n_cov);
  double **psiPro = doubleMatrix(n_covp, n_covp);
  double **psi1Pro = doubleMatrix(n_covp, n_covp);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      deltaPro[i][j] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < n_covp; j++)
    for (i = 0; i < n_covp; i++)
      psiPro[i][j] = psiVar[itemp++];

  if (*unconst == 1) {
    itemp = 0;
    for (j = 0; j < n_covp; j++)
      for (i = 0; i < n_covp; i++)
	psi1Pro[i][j] = psiVar[itemp++];
  }

  /** known Z star **/
  for (i = 0; i < *n_samp; i++) {
    if ((treat[i] == 1) && (Y[i] == (*J+1))) {
      Zstar[i] = 1;
      Y0[i] = *J;
    } else if ((treat[i] == 1) && (Y[i] == 0)) {
      Zstar[i] = 0;
      Y0[i] = 0;
    } else { /* random draw if not known */
      Zstar[i] = (unif_rand() < 0.5);
    }
    if (*unconst == 2) 
      X[i][n_covp-1] = Zstar[i];
  }
  
  /** MCMC **/
  itemp = 0; deltaCounter[0] = 0; psiCounter[0] = 0; psi1Counter[0] = 0;
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... *** \n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    itempZstar0 = 0; itempZstar1 = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar */
      if ((treat[i] == 0) || ((treat[i] == 1) && (Y[i] < (*J + 1)) && (Y[i] > 0))) {
	Xdelta[i] = 0;  Xpsi[i] = 0;  Xpsi1[i] = 0;
	for (j = 0; j < *n_cov; j++) { 
	  Xdelta[i] += X[i][j]*delta[j];
	  Xpsi[i] += X[i][j]*psi[j];
	  if (*unconst == 1) {
	    Xpsi1[i] += X[i][j]*psi1[j];
	  } 
	}
	if (*unconst == 2) {
	  Xpsi1[i] = Xpsi[i] + psi[n_covp-1];
	}
	if (*unconst > 0) {
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
	if (*unconst == 1) {
	  if (Zstar[i] == 1) {
	    Y0temp0[itempZstar0] = Y0[i];
	    for (j = 0; j < n_covp; j++) {
	      Xtemp0[itempZstar0][j] = X[i][j];
	    }
	    itempZstar0++;
	  } else {
	    Y0temp1[itempZstar1] = Y0[i];
	    for (j = 0; j < n_covp; j++) {
	      Xtemp1[itempZstar1][j] = X[i][j];
	    }
	    itempZstar1++;
	  }	    
	}
	if (*unconst == 2) { 
	  X[i][n_covp-1] = Zstar[i];
	}
      }
    }
    
    /* Sample delta */
    BinomLogit(Zstar, X, delta, *n_samp, 1, *n_cov, delta0, A0delta, deltaPro, 1, deltaCounter);
    
    /* Sample psi */
    if (*unconst == 1) {
      BinomLogit(Y0temp0, Xtemp0, psi, itempZstar0, *J, n_covp, psi0, A0psi, psiPro, 1, psiCounter);
      BinomLogit(Y0temp1, Xtemp1, psi1, itempZstar1, *J, n_covp, psi01, A0psi1, psi1Pro, 1, psi1Counter);
    } else {
      BinomLogit(Y0, X, psi, *n_samp, *J, n_covp, psi0, A0psi, psiPro, 1, psiCounter);
    }

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = delta[j];
	}      
	for (j = 0; j < n_covp; j++) {
	  allresults[itemp++] = psi[j];
	}
	if (*unconst == 1) {
	  for (j = 0; j < n_covp; j++) {
	    allresults[itemp++] = psi1[j];
	  }
	}
	allresults[itemp++] = ((double) *deltaCounter / (double) (main_loop + 1));
	allresults[itemp++] = ((double) *psiCounter / (double) (main_loop + 1));
	if (*unconst == 1) {
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
	Rprintf("    Metropolis acceptance ratio (sensitive): %3g\n", 
		fprec((double) *deltaCounter / (double) (main_loop + 1), 3));
	if (*unconst == 1) {
	  Rprintf("    Metropolis acceptance ratios (control): %3g (psi)  %3g (psi1)\n", 
		  fprec((double) *psiCounter / (double) (main_loop + 1), 3),
		  fprec((double) *psi1Counter / (double) (main_loop + 1), 3));
	} else {
	  Rprintf("    Metropolis acceptance ratio (control): %3g\n", 
		  fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	}
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
  FreeMatrix(A0psi, n_covp);
  FreeMatrix(A0psi1, n_covp);
  FreeMatrix(deltaPro, *n_cov);
  FreeMatrix(psiPro, n_covp);
  FreeMatrix(psi1Pro, n_covp);
}


/** 
   Item Count Technique Binomial Regression for the Multiple Sensitive Item Design
**/


void ictregBinomMulti(int *Y,             /* outcome vector */
		      int *J,             /* # of control items */
		      int *n_samp,        /* sample size */
		      int *n_draws,       /* # of MCMC draws */
		      int *treat,         /* treatment indicator vector: 0, ..., tmax */
		      int *tmax,          /* number of sensitive items */
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

  int i, j, k, main_loop, itemp, itempS, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int **deltaCounter = intMatrix(*tmax, 1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);        /* acceptance ratio for psi */
  double dtemp1, dtemp2;

  /** get random seed **/
  GetRNGstate();

  /** Parameters for sensitive items **/
  int n_dim = *n_cov;
  if (*unconst) { /* dimension of delta */
    n_dim = *n_cov + 1;
  }

  double **deltaMatrix = doubleMatrix(*tmax, n_dim);

  itemp = 0; 
  for (i = 0; i < *tmax; i++) 
    for (j = 0; j < n_dim; j++)
      deltaMatrix[i][j] = delta[itemp++];

  /* PdoubleMatrix(deltaMatrix, *tmax, n_dim);
     R_FlushConsole(); */

  /** Data **/
  int *Zstar = intArray(*n_samp);
  int *Y0 = intArray(*n_samp);
  int *Zstartemp = intArray(*n_samp);
  double *Xdelta = doubleArray(*n_samp);
  double *Xdelta1 = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim); 
  double **Xtemp = doubleMatrix(*n_samp, n_dim); 
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  for (i = 0; i < *n_samp; i++)
    Y0[i] = Y[i];

  /** Prior **/
  double ***A0delta = doubleMatrix3D(*tmax, n_dim, n_dim);
  double **A0psi = doubleMatrix(*n_cov, *n_cov);
  double **m0delta = doubleMatrix(*tmax, n_dim);

  itemp = 0;
  for (i = 0; i < *tmax; i++)
    for (k = 0; k < n_dim; k++)
      for (j = 0; j < n_dim; j++)
	A0delta[i][j][k] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  itemp = 0;
  for (i = 0; i < *tmax; i++)
    for (j = 0; j < n_dim; j++)
      m0delta[i][j] = delta0[itemp++];

  /* PdoubleMatrix3D(A0delta, *tmax, n_dim, n_dim);
     PdoubleMatrix(m0delta, *tmax, n_dim); 
     PdoubleMatrix(A0psi, *n_cov, *n_cov);
     R_FlushConsole(); */

  /** Proposal precisoin **/
  double ***deltaPro = doubleMatrix3D(*tmax, n_dim, n_dim);
  double **psiPro = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (i = 0; i < *tmax; i++)
    for (k = 0; k < n_dim; k++)
      for (j = 0; j < n_dim; j++)
	deltaPro[i][j][k] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      psiPro[i][j] = psiVar[itemp++];

  /* PdoubleMatrix3D(deltaPro, *tmax, n_dim, n_dim); 
     PdoubleMatrix(psiPro, *n_cov, *n_cov);
     R_FlushConsole(); */

  /** MCMC **/
  itempS = 0; psiCounter[0] = 0; 
  for (i = 0; i < *tmax; i++) {
    deltaCounter[i][0] = 0; 
  }
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... *** \n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar for treated units */
      if (treat[i] > 0) {
	if (Y[i] == (*J+1)) {
	  Zstar[i] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[i] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0;  Xpsi[i] = 0;  
	  for (j = 0; j < *n_cov; j++) {
	    Xdelta[i] += X[i][j] * deltaMatrix[treat[i]-1][j];
	    Xpsi[i] += X[i][j] * psi[j];
	  }
	  if (*unconst) {
	    Xdelta1[i] = Xdelta[i] + (Y[i] - 1) * deltaMatrix[treat[i]-1][*n_cov];
	    Xdelta[i] += Y[i] * deltaMatrix[treat[i]-1][*n_cov];
	    dtemp1 = exp(Xdelta1[i] - log1p(exp(Xdelta1[i])) + 
			 dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  } else {
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  } 
	  if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	    Zstar[i] = 1;
	    Y0[i] = Y[i] - 1;
	  } else { 
	    Zstar[i] = 0;
	    Y0[i] = Y[i];
	  }
	}
	if (*unconst) {
	  X[i][*n_cov] = Y0[i];
	}
      }
    }

    /* Sample delta */
    for (k = 1; k <= *tmax; k++) {
      itemp = 0;
      for (i = 0; i < *n_samp; i++) {
	if (treat[i] == k) {
	  Zstartemp[itemp] = Zstar[i];
	  for (j = 0; j < n_dim; j++) {
	    Xtemp[itemp][j] = X[i][j];
	  }
	  itemp++;
	}
      }
      BinomLogit(Zstartemp, Xtemp, deltaMatrix[k-1], itemp, 1, n_dim, 
		 m0delta[k-1], A0delta[k-1], deltaPro[k-1], 1, deltaCounter[k-1]);
    }

    /* Sample psi */
    BinomLogit(Y0, X, psi, *n_samp, *J, *n_cov, psi0, A0psi, psiPro, 1, psiCounter);

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	for (k = 0; k < *tmax; k++) {
	  for (j = 0; j < n_dim; j++) {
	    allresults[itempS++] = deltaMatrix[k][j];
	  }      
	}
	for (j = 0; j < *n_cov; j++) {
	  allresults[itempS++] = psi[j];
	}
	for (k = 0; k < *tmax; k++) {
	  allresults[itempS++] = ((double) deltaCounter[k][0] / (double) (main_loop + 1));
	}
	allresults[itempS++] = ((double) *psiCounter / (double) (main_loop + 1));
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
	Rprintf("    Metropolis acceptance ratio (sensitive): %3g", 
		fprec((double) deltaCounter[0][0] / (double) (main_loop + 1), 3));
	for (k = 1; k < *tmax; k++) {
	  Rprintf(" %3g", fprec((double) deltaCounter[k][0] / (double) (main_loop + 1), 3));
	}
	Rprintf("\n    Metropolis acceptance ratio (control): %3g\n", 
		fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(psiCounter);
  free(Zstar);
  free(Y0);
  free(Zstartemp);
  free(Xdelta);
  free(Xdelta1);
  free(Xpsi);
  free(m0delta);
  FreeintMatrix(deltaCounter, *tmax);
  FreeMatrix(deltaMatrix, *tmax);
  FreeMatrix(X, *n_samp);
  FreeMatrix(Xtemp, *n_samp);
  Free3DMatrix(A0delta, *tmax, n_dim);
  FreeMatrix(A0psi, *n_cov);
  Free3DMatrix(deltaPro, *tmax, n_dim);
  FreeMatrix(psiPro, *n_cov);
}


/** 
  Item Count Technique Binomial Mixed Effects Regression for the Standard Design
  
  for now, constrained model only

**/

void ictregBinomMixed(int *Y,             /* outcome vector */
		      int *J,             /* # of control items */
		      int *n_samp,        /* sample size */
		      int *n_draws,       /* # of MCMC draws */
		      int *treat,         /* treatment indicator vector: 0 or 1 */
		      double *Xall,       /* fixed effects covariates in a vector form */
		      double *delta,      /* fixed effects coefs for sensitive item */
		      double *psi,        /* fixed effects coefs for control items */ 
		      int *n_cov,         /* # of fixed effects covariates */
		      double *delta0,     /* prior mean for delta */
		      double *psi0,       /* prior mean for psi */
		      double *A0deltaAll, /* prior precision for delta */
		      double *A0psiAll,   /* prior precision for psi */
		      double *deltaVar,   /* proposal variance for delta (fixed effects) */
		      double *psiVar,     /* proposal variance for psi (fixed effects) */
		      int *grp,           /* group indicator, 0, 1, ..., G-1 */
		      int *n_grp,         /* number of groups, G */
		      int *max_samp_grp,  /* max # of obs within each group */
		      double *Xall_rand,  /* random effects covariates */
		      int *n_rand,        /* # of random effects covariates */
		      double *tune_gamma, /* tuning constants for random effects (sensitive) */
		      double *tune_zeta,  /* tuning constants for random effects (control) */
		      double *SigmaAll,   /* covariance for random effects (sensitive) */
		      double *PhiAll,     /* covariance for random effects (control) */
		      int *s0,            /* prior df for Sigma (sensitive) */
		      double *S0All,      /* prior scale matrix for Sigma (sensitive) */
		      int *t0,            /* prior df for Phi (control) */
		      double *T0All,      /* prior scale matrix for Phi (control) */
		      int *burnin,        /* number of burnins */
		      int *keep,          /* keep every *th draw */
		      int *verbose,       /* want to print progress? */
		      double *allresults  /* storage for all results */
		      ) {

  int i, j, k, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int *deltaCounter = intArray(1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);    /* acceptance ratio for psi */
  int *gammaCounter = intArray(*n_grp);
  int *zetaCounter = intArray(*n_grp);
  int *vitemp = intArray(*n_grp);
  double dtemp1, dtemp2;

  /** get random seed **/
  GetRNGstate();

  /** parameters for random effects **/
  double **gamma = doubleMatrix(*n_grp, *n_rand); /* sensitive */
  double *gamma0 = doubleArray(*n_rand);
  double **Sigma = doubleMatrix(*n_rand, *n_rand);
  double **zeta = doubleMatrix(*n_grp, *n_rand);  /* control */
  double *zeta0 = doubleArray(*n_rand);
  double **Phi = doubleMatrix(*n_rand, *n_rand);

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      Sigma[i][j] = SigmaAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      Phi[i][j] = PhiAll[itemp++];

  /* starting values for random effects from prior */
  for (j = 0; j < *n_rand; j++) {
    gamma0[j] = 0; 
    zeta0[j] = 0;
  }
  for (j = 0; j < *n_grp; j++) { 
    rMVN(gamma[j], gamma0, Sigma, *n_rand);
    rMVN(zeta[j], zeta0, Phi, *n_rand);
  }
  
  /** Data **/
  int *Zstar = intArray(*n_samp);
  int *Y0 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, *n_cov); 
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double ***Vgrp = doubleMatrix3D(*n_grp, *max_samp_grp, *n_rand);
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++) {
    vitemp[j] = 0;
  }
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_rand; j++)
      Vgrp[grp[i]][vitemp[grp[i]]][j] = Xall_rand[itemp++];
    vitemp[grp[i]]++;
  }
 
  /* PdoubleMatrix3D(Vgrp, *n_grp, *max_samp_grp, *n_rand); */

  /** Prior **/
  double **A0delta = doubleMatrix(*n_cov, *n_cov);
  double **A0psi = doubleMatrix(*n_cov, *n_cov);
  double **S0 = doubleMatrix(*n_rand, *n_rand);
  double **T0 = doubleMatrix(*n_rand, *n_rand);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0delta[i][j] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      S0[i][j] = S0All[itemp++];

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      T0[i][j] = T0All[itemp++];

  /** Proposal precisoin **/
  double **deltaPro = doubleMatrix(*n_cov, *n_cov);
  double **psiPro = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      deltaPro[i][j] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      psiPro[i][j] = psiVar[itemp++];

  /** MCMC **/
  itemp = 0; deltaCounter[0] = 0; psiCounter[0] = 0; 
  for (j = 0; j < *n_grp; j++) {
    gammaCounter[j] = 0;
    zetaCounter[j] = 0;
  }
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... ***\n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar */
      if ((treat[i] == 1) && (Y[i] == (*J+1))) {
	Zstar[i] = 1;
	Y0[i] = *J;
      } else if ((treat[i] == 1) && (Y[i] == 0)) {
	Zstar[i] = 0;
	Y0[i] = 0;
      } else {
	Xdelta[i] = 0;  Xpsi[i] = 0;  
	for (j = 0; j < *n_cov; j++) { /* fixed effects */
	  Xdelta[i] += X[i][j]*delta[j];
	  Xpsi[i] += X[i][j]*psi[j];
	}
	for (j = 0; j < *n_rand; j++) { /* random effects */
	  Xdelta[i] += Vgrp[grp[i]][vitemp[grp[i]]][j] * gamma[grp[i]][j];
	  Xpsi[i] += Vgrp[grp[i]][vitemp[grp[i]]][j] * zeta[grp[i]][j];
	}
	dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
		     dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	  Zstar[i] = 1;
	  Y0[i] = Y[i] - treat[i];
	} else { 
	  Zstar[i] = 0;
	  Y0[i] = Y[i];
	}
      }
      vitemp[grp[i]]++;
    }

    /* Sample delta */
    BinomLogitMixed(Zstar, X, Vgrp, grp, delta, gamma, Sigma, *n_samp, 1, *n_cov, *n_rand, *n_grp, 
		    delta0, A0delta, *s0, S0, deltaPro, tune_gamma, 1, deltaCounter, gammaCounter);
    
    /* Sample psi */
    BinomLogitMixed(Y0, X, Vgrp, grp, psi, zeta, Phi, *n_samp, *J, *n_cov, *n_rand, *n_grp, 
		    psi0, A0psi, *t0, T0, psiPro, tune_zeta, 1, psiCounter, zetaCounter);

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	/* fixed effects */
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = delta[j];
	}      
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = psi[j];
	}
	/* Var-Cov matrices */
	for (j = 0; j < *n_rand; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = Sigma[j][k];
	  }      
	}
	for (j = 0; j < *n_rand; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = Phi[j][k];
	  }      
	}
	/* random effects */
	for (j = 0; j < *n_grp; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = gamma[j][k];
	  }
	}
	for (j = 0; j < *n_grp; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = zeta[j][k];
	  }
	}
	/* acceptance ratios */
	allresults[itemp++] = ((double) *deltaCounter / (double) (main_loop + 1));
	allresults[itemp++] = ((double) *psiCounter / (double) (main_loop + 1));
	for (j = 0; j < *n_rand; j++)
	  allresults[itemp++] = ((double) gammaCounter[j] / (double) (main_loop + 1));
	for (j = 0; j < *n_rand; j++)
	  allresults[itemp++] = ((double) zetaCounter[j] / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("\n    Metropolis acceptance ratio (fixed effects): %3g (sensitive) %3g (control)", 
		fprec((double) *deltaCounter / (double) (main_loop + 1), 3),
		fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	Rprintf("\n    Metropolis acceptance ratio (random effects: sensitive):");
	for (j = 0; j < *n_grp; j++) 
	  Rprintf(" %3g", fprec((double) gammaCounter[j] / (double) (main_loop + 1), 3));
	Rprintf("\n    Metropolis acceptance ratio (random effects: control):");
	for (j = 0; j < *n_grp; j++) 
	  Rprintf(" %3g", fprec((double) zetaCounter[j] / (double) (main_loop + 1), 3));
	Rprintf("\n");
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
  free(Xdelta);
  free(Xpsi);
  free(deltaCounter);
  free(psiCounter);
  free(gammaCounter);
  free(zetaCounter);
  free(vitemp);
  free(gamma0);
  free(zeta0);
  FreeMatrix(Sigma, *n_rand);
  FreeMatrix(Phi, *n_rand);
  FreeMatrix(X, *n_samp);
  FreeMatrix(A0delta, *n_cov);
  FreeMatrix(A0psi, *n_cov);
  FreeMatrix(S0, *n_rand);
  FreeMatrix(T0, *n_rand);
  FreeMatrix(deltaPro, *n_cov);
  FreeMatrix(psiPro, *n_cov);
  Free3DMatrix(Vgrp, *n_grp, *max_samp_grp);
}



/**
   A Random Walk Metroplis Sampler for Binomial Logistic Regression 
   with Independent Normal Prior
   
   proposal distribution is the univariate normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.
**/

void BinomLogit(int *Y,        /* outcome variable: 0, 1, ..., J */
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
} /* end of BinomLogit */


/**
   A Random Walk Metroplis Sampler for Binomial Logistic Mixed Effects 
   Regression with Independent Normal Prior and Normal random effects.
   
   proposal distribution for fixed effects is the normal whose mean is
   the current value and variance is given by the input. 

   proposal distribution for random effects is the multivariate normal
   whose mean is the current value and variance is given by the
   current value of covariance matrix multiplied by the input tuning
   parameter. 

**/

void BinomLogitMixed(int *Y,          /* outcome variable: 0, 1, ..., J */
		     double **X,      /* (N x K) covariate matrix for
					 fixed effects */
		     double ***Z,     /* covariates for random effects 
				         organized by groups */
		     int *grp,        /* group indicator, 0, 1, ..., G-1 */
		     double *beta,    /* K coefficients for fixed effects */
		     double **gamma,  /* (G x L) matrix of random effects */
		     double **Psi,    /* LxL precision matrix for random effecs */
		     int n_samp,      /* # of obs */
		     int J,           /* size of binomial, J */
		     int n_fixed,     /* # of fixed effects, K */
		     int n_random,    /* # of random effects, L */
		     int n_grp,       /* # of groups, G */
		     double *beta0,   /* K dimensional prior mean vector */
		     double **A0,     /* (K x K) prior precision */
		     int tau0,        /* prior df for Psi */
		     double **T0,     /* prior scale for Psi */
		     double **tune_fixed,  /* K proposal variance-covariance matrix */
		     double *tune_random, /* tuning constant for random effects of each group */
		     int n_gen,        /* # of MCMC draws */
		     int *acc_fixed,   /* # of acceptance for fixed effects */
		     int *acc_random   /* # of acceptance for random effects */
		     ) {
  
  int i, j, k, main_loop;
  int *vitemp = intArray(n_grp);
  double numer, denom;
  /* proposal values */
  double *beta1 = doubleArray(n_fixed);
  double *gamma1 = doubleArray(n_random);
  /* prior for gamma = 0 */
  double *gamma0 = doubleArray(n_random);
  /* data holders */
  double *Xbeta = doubleArray(n_samp);
  double *Xbeta1 = doubleArray(n_samp);
  double *Zgamma = doubleArray(n_samp);
  double *Zgamma1 = doubleArray(n_samp);
  /* matrix holders */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  for (j = 0; j < n_fixed; j++)
    beta1[j] = beta[j];

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /** initializing Xbeta and Zgamma **/
  for (j = 0 ; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    Xbeta[i] = 0; Zgamma[i] = 0;
    for (j = 0; j < n_fixed; j++) { 
      Xbeta[i] += X[i][j] * beta[j];
    }
    Xbeta1[i] = Xbeta[i];
    for (j = 0; j < n_random; j++) {
      Zgamma[i] += Z[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
    }
    Zgamma1[i] = Zgamma[i];
    vitemp[grp[i]]++;
  }

  /** MCMC Sampler starts here **/
  for (main_loop = 0; main_loop < n_gen; main_loop++) {

    /** STEP 1: Update Fixed Effects **/
    rMVN(beta1, beta, tune_fixed, n_fixed);
    /* Calculating the ratio (log scale) */
    /* prior */
    numer = dMVN(beta1, beta0, A0, n_fixed, 1);
    denom = dMVN(beta, beta0, A0, n_fixed, 1);   
    /* likelihood */
    for (i = 0; i < n_samp; i++) {
      Xbeta1[i] = Xbeta[i];
      for (j = 0; j < n_fixed; j++) {
	Xbeta1[i] -= X[i][j] * (beta[j] - beta1[j]);
      }
      denom += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma[i])), 1);
      numer += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta1[i]-Zgamma[i])), 1);
    }
    /* Rejection */
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      acc_fixed[0]++;
      for (j = 0; j < n_fixed; j++)
	beta[j] = beta1[j];
      for (i = 0; i < n_samp; i++) 
	Xbeta[i] = Xbeta1[i];
    }
    /** STEP 2: Update Random Effects Given Fixed Effects **/
    dinv(Psi, n_random, mtemp);
    for (i = 0; i < n_random; i++)
      for (j = 0; j < n_random; j++)
	mtemp[i][j] *= tune_random[j];
    for (j = 0; j < n_grp; j++) {
      rMVN(gamma1, gamma[j], mtemp, n_random);
      /* Calculating the ratio (log scale) */
      /* prior */
      numer = dMVN(gamma1, gamma0, Psi, n_random, 1);
      denom = dMVN(gamma[j], gamma0, Psi, n_random, 1); 
      /* likelihood for group j */
      for (k = 0 ; k < n_grp; k++)
	vitemp[k] = 0;
      for (i = 0; i < n_samp; i++) {
	if (grp[i] == j) {
	  for (k = 0; k < n_random; k++)
	    Zgamma1[i] = Zgamma[i] - Z[j][vitemp[j]][k]*(gamma[j][k]-gamma1[k]);
	  denom += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma[i])), 1);
	  numer += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma1[i])), 1);
	}
	vitemp[grp[i]]++;
      }
      /* Rejection */
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	acc_random[j]++;
	for (k = 0; k < n_random; k++)
	  gamma[j][k] = gamma1[k];
	for (i = 0; i < n_samp; i++) {
	  if (grp[i] == j) {
	    Zgamma[i] = Zgamma1[i];
	  }      
	}
      }
    }
    
    /** STEP 3: Update Psi **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (i = 0; i < n_grp; i++)
      for (j = 0; j < n_random; j++)
	for (k = 0; k < n_random; k++)
	  mtemp[j][k] += gamma[i][j] * gamma[i][k];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, tau0+n_grp, n_random);
  }

  /* freeing memory */
  free(beta1);
  free(gamma1);
  free(gamma0);
  free(vitemp);
  free(Xbeta);
  free(Xbeta1);
  free(Zgamma);
  free(Zgamma1);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects logit */
