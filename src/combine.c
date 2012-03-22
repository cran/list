#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

/** 
   Combining Endorsement Experiments with List Experiments (Multiple
   Sensitive Items) 

   Constrained model with Robit regression for the sensitive item

**/



/**
  
  ORIGINAL: endorsement experiment model code from endorse package
  where bNormalReg needs to be changed to bNormalReg2

 **/

void CombineEndorseList(/* Starging the endorsement code stuff */
		int *dY,          /* ordinal outcome variable: 0, 1,
				     ..., J-1. Length N * J */
		int *dT,          /* endorsement matrix. Length N * J */
		double *dZ,      /* covariate matirx, length N * M */
		/*   Data Structure   */
		int *n_samp,     /* # of obs: N */ 
		int *n_pol,     /* # of policies: J */
		int *n_dim,     /* dimension of covariates: M  */
		int *n_act,      /* # of actors: K */
		int *n_cat,      /* # of categories for each questions:
				    L_j. Length  J*/
		int *max_n_cat,  /* max # of categories */
		/*   Starting values   */
		double *X,      /* vector of ideal points. 
				    No need to include constant
				    Length N */
		double *dS,       /* matrix of support, s_{ij}(k)
				     Length (N * J) */
		double *sbeta,    /* (alpha_j, beta_j), length J * 2 */
		double *stau,     /* cut points for each question;
				     the first cut point is set to 0
				     and the last one set to tau_{L_j - 1}+1000:
				     Length J * max {L_j} */
		double *slambda,  /* lambda, length (J * M) * K  */
		double *somega2,   /* omega2: variance of s given lambda,
				     length J * K */
		double *dtheta,   /* theta, vector of length K * M  */
		double *phi2,     /* phi, vector of length K * M */
		double *delta,    /* delta, vector of length M  */
		/*  Prior parameters  */
		double *dmu_beta,   /* prior mean of factor loadings */
		double *dA0_beta,     /* prior precision of factor loadings,
				    length 2 * 2: can be set to zero to
				    induce improper prior for beta alone */
		double *dA0_x,     /* known prior precision of ideal points, identical for all i */
		double *dmu_theta, /* prior mean for theta, length K */
		double *dA0_theta, /* prior precision of theta, identical for all k */
		double *mu_delta,  /* prior mean of delta, length M */
		double *dA0_delta,  /* prior precision of delta, length M * M */
		double *s0_omega2,      /* prior scale for InvChi2 of omega */
		int *nu0_omega2,         /* prior d.f. for InvChi2 of omega */
		double *s0_phi2,    /* prior scale for InvChi2 of phi */
		int *nu0_phi2,      /* prior d.f. for InvChi2 of phi */
		/* MCMC settings */
		int *n_gen,       /* # of gibbs draws */
		int *burn,       /* # of burnin period */
		int *thin,       /* thinning */
		int *mda,        /* marginal data augmentation? */
		int *mh,         /* Metropolis-Hasting step? */
		double *prop,    /* proposal variance for MH step */
		int *x_sd,       /* Output of x: sd if 1, sample if 0 */
		int *tau_out,
		double *betaStore,
		double *tauStore,
		double *xStore,
		double *lambdaStore,
		double *thetaStore,
		double *deltaStore,
		double *accept_ratio,
		/* list experiment stuff */
		int *Ylist,         /* outcome vector */
		int *J,             /* # of control items */
		int *treat,         /* treatment indicator vector: 0, ..., tmax */
		double *psi,        /* coefs for control items */ 
		double *psi0,       /* prior mean for psi */
		double *A0psiAll,   /* prior precision for psi */
		double *psiVar,     /* proposal variance for psi */
		int *df,             /* degrees of freedom for t-distribution */
		int *ceiling,       /* ceiling effects */
		int *floor,         /* floor effects */
		int *verbose,       /* want to print progress? */
		double *allresults  /* storage for results */
	       ){
  /* loop counters */
  int i, j, k, l, m, n, main_loop, itemp;
  int ibeta = 0, itau = 0, ix = 0, ilambda = 0, itheta = 0, idelta = 0, keep = 1;
  double varx;

  /* setting up the endorsement stuff */
  /* storage vectors */
  int *Y_j = intArray(*n_samp);
  double *beta = doubleArray(2);
  double *mu_beta = doubleArray(2);
  double *s = doubleArray(1);
  double *var_epsilon = doubleArray(1);
  var_epsilon[0] = 1;
  double *mu_s = doubleArray(1);
  double *x_i = doubleArray(1);
  double *mu_x = doubleArray(1);
  double *lambda_jk = doubleArray(*n_dim);
  double *omega2_jk = doubleArray(1);
  double *mu_lambda = doubleArray(*n_dim);
  double *stemp = doubleArray(*n_samp);
  double *ztemp = doubleArray(*n_samp * *n_dim);
  double *theta_km = doubleArray(*n_dim);
  double *mu_theta = doubleArray(1);
  double *sig2_x = doubleArray(1);
  int *accept = intArray(*n_pol);
  for (j = 0; j < *n_pol; j++)
    accept[j] = 0;
  int *temp_accept = intArray(1);

  /* storage matrices */
  /** data **/
  int **T = intMatrix(*n_samp, *n_pol); /* treatment */
  int **Y = intMatrix(*n_samp, *n_pol); /* response */
  double **Z = doubleMatrix(*n_samp, *n_dim); /* covariates */

  /** starting values **/
  double **S = doubleMatrix(*n_samp, *n_pol); /* support parameter */
  double **Beta = doubleMatrix(*n_pol, 2); /* alpha and beta */
  double **Tau = doubleMatrix(*n_pol, *max_n_cat); /* cut point */
  double ***Lambda = doubleMatrix3D(*n_pol, *n_act, *n_dim); /* lambda */
  double **Omega2 = doubleMatrix(*n_pol, *n_act); /* omega2 */
  double **Theta = doubleMatrix(*n_act, *n_dim); /* theta */

  /** prior mean and precision **/
  double **Mu_beta = doubleMatrix(*n_pol, 2);
  double **A0_beta = doubleMatrix(2, 2);
  double **A0_s = doubleMatrix(1, 1);
  double **A0_x = doubleMatrix(1, 1);
  double **A0_lambda = doubleMatrix(*n_dim, *n_dim);
  double **A0_theta = doubleMatrix(*n_dim, *n_dim);
  double **A0_delta = doubleMatrix(*n_dim, *n_dim);

  /** matrices for regression **/
  double **Ystar = doubleMatrix(*n_samp, *n_pol); /* latent random utility */
  double **U = doubleMatrix(*n_samp+2, 3); /* x_{i} + s_{ij}(k) */
  double **D_s = doubleMatrix(2, 2);
  double **D_x = doubleMatrix(*n_pol+1, 2);
  double **D_theta = doubleMatrix(*n_pol+1, 2);
  double **D_delta = doubleMatrix(*n_samp + *n_dim, *n_dim + 1);

  /* get random seed */
  GetRNGstate();

  /* packing data */
  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for (i = 0; i < *n_samp; i++)
      Y[i][j] = dY[itemp++];

  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for(i = 0; i < *n_samp; i++)
      T[i][j] = dT[itemp++];

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (i = 0; i < *n_samp; i++)
      Z[i][m] = dZ[itemp++];

  /* packing starting values */
  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for (i = 0; i < *n_samp; i++)
      S[i][j] = dS[itemp++];

  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < *n_pol; j++)
      Beta[j][m] = sbeta[itemp++];

  itemp = 0;
  for (l = 0; l < *max_n_cat; l++)
    for (j = 0; j < *n_pol; j++)
      Tau[j][l] = stau[itemp++];

  itemp = 0;
  for (k = 0; k < *n_act; k++)
    for (j = 0; j < *n_pol; j++)
      for (m = 0; m < *n_dim; m++)
      Lambda[j][k][m] = slambda[itemp++];

  itemp = 0;
  for (k = 0; k < *n_act; k++)
    for (j = 0; j < *n_pol; j++)
      Omega2[j][k] = somega2[itemp++];

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (k = 0; k < *n_act; k++)
      Theta[k][m] = dtheta[itemp++];

  /* packing prior mean */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < *n_pol; j++)
      Mu_beta[j][m] = dmu_beta[itemp++];

  /* packing prior precision */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (n = 0; n < 2; n++)
      A0_beta[n][m] = dA0_beta[itemp++];

  for (m = 0; m < *n_dim; m++)
    for (n = 0; n < *n_dim; n++)
      A0_lambda[m][n] = 0;

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (n = 0; n < *n_dim; n++)
      A0_theta[m][n] = dA0_theta[itemp++];

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (n = 0; n < *n_dim; n++)
      A0_delta[n][m] = dA0_delta[itemp++];

  A0_x[0][0] = *dA0_x;
  sig2_x[0] = (1 / *dA0_x);


  /** Setting up the list experiment stuff **/
  int progress = 1;
  int itempS;
  int *treatSum = intArray(*n_act);
  int *psiCounter = intArray(1);        /* acceptance ratio for psi */
  int **Zstar = intMatrix(*n_act, *n_samp);
  int *Y0 = intArray(*n_samp);
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double ***Xtemp = doubleMatrix3D(*n_act, *n_samp + *n_pol * *n_dim, *n_dim + 1); 
  double **Prec = doubleMatrix(*n_act, *n_samp + *n_pol * *n_dim); 
  double **A0psi = doubleMatrix(*n_dim, *n_dim);
  double dtemp, dtemp1, dtemp2, dtemp3;
  double **psiPro = doubleMatrix(*n_dim, *n_dim);

  itemp = 0;
  for (j = 0; j < *n_dim; j++)
    for (i = 0; i < *n_dim; i++)
      psiPro[i][j] = psiVar[itemp++];

  for (i = 0; i < *n_act; i++)
    treatSum[i] = 0;

  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Ylist[i];
    if (treat[i] > 0) {
      for (j = 0; j < *n_dim; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = Z[i][j];
      treatSum[treat[i]-1]++;
      if ((ceiling[treat[i]-1] == 1) && (Ylist[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Ylist[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
  }

  itemp = 0;
  for (j = 0; j < *n_dim; j++)
    for (i = 0; i < *n_dim; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  /** Gibbs Sampler **/
  itempS;
  Rprintf("Start Gibbs Sampler\n");
  for (main_loop = 1; main_loop <= *n_gen; main_loop++) {
    Rprintf("%5d\n", main_loop);
    if (*verbose)
      if (main_loop % 100 == 0)
	Rprintf("%6d / %6d\n", main_loop, *n_gen);

    /** Endorsement experiment model sampling begins here **/
    /** start sampling alpha and beta **/
    for (j = 0; j < *n_pol; j++) {

      /*** vectors ***/
      double *tau = doubleArray(n_cat[j]);
      double *MHprop = doubleArray(n_cat[j]-2);

      /*** proposal variance vector ***/
      for (l = 0; l < (n_cat[j]-2); l++)
	MHprop[l] = *prop;

      /*** response of each question ***/
      for (i = 0; i < *n_samp; i++)
	Y_j[i] = Y[i][j];
      
      /*** systematic component of utility ***/
      for (i = 0; i < *n_samp; i++) {
	U[i][0] = -1;
	U[i][1] = X[i] + S[i][j];
      }

      /*** starting values of alpha and beta ***/
      for (m = 0; m < 2; m++)
	beta[m] = Beta[j][m];

      /*** starting values of tau ***/
      for (l = 0; l < n_cat[j]; l++)
	tau[l] = Tau[j][l];

      
      /*** prior mean ***/
      for (m = 0; m < 2; m++)
	mu_beta[m] = Mu_beta[j][m];

      /*** set acceptance indicator to 0 ***/
      temp_accept[0] = 0;

      /*** draw posterior ***/
      endorseoprobitMCMC(Y_j, U, beta, tau, *n_samp, 2, n_cat[j], 1,
			  mu_beta, A0_beta, *mda, *mh, MHprop, temp_accept, 1);

      /*** update acceptance counter ***/
      accept[j] += temp_accept[0];
      accept_ratio[j] = (double) accept[j] / (double) main_loop;

      /*** packing sampled alpha and beta ***/
      for (m = 0; m < 2; m++)
	Beta[j][m] = beta[m];

      /*** packing sampled tau ***/
      for (l = 0; l < n_cat[j]; l++)
	Tau[j][l] = tau[l];
      
      /*** packing sampled latent random utility ***/
      for (i = 0; i < *n_samp; i++)
	Ystar[i][j] = U[i][2];

      /*** storing alpha and beta ***/
      if(main_loop > *burn) {
	if(keep == *thin) {

	  for (m = 0; m < 2; m++)
	    betaStore[ibeta++] = Beta[j][m];

	  if (*tau_out) {
	    for (l = 0; l < (*max_n_cat-1); l++)
	      tauStore[itau++] = Tau[j][l];
	  }
	}
      }
      
      /** print and store acceptance ratio  **/      
      if (*mh) {
	if (main_loop % 100 == 0) {
	  Rprintf("Acceptance rate for Tau.%2d: %4g\n", (j + 1), accept_ratio[j]);
	}
      }

      R_FlushConsole();
      R_CheckUserInterrupt();
    } /** end of sampling alpha and beta **/

    /** start sampling s **/
    for (i = 0; i < *n_samp; i++){
      for (j = 0; j < *n_pol; j++){

	k = T[i][j];
	
	/*** if not control, sample s ***/
	if (k > 0) {

	  D_s[0][0] = Beta[j][1];
	  D_s[0][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*X[i];

	  mu_s[0] = 0;
	  for (m = 0; m < *n_dim; m++)
	    mu_s[0] += Lambda[j][k-1][m] * Z[i][m];

	  A0_s[0][0] = (1 / Omega2[j][k-1]);

	  bNormalReg2(D_s, s, var_epsilon, 1, 1, 1, 1, mu_s, A0_s,
		     1, 1, 1, 1);

	  /*** packing sampled s ***/
	  S[i][j] = s[0];
	}

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    }/** end of sampling s **/

    /** start sampling x **/
    for (i = 0; i < *n_samp; i++) {

      for (j = 0; j < *n_pol; j++) {
	D_x[j][0] = Beta[j][1];
	D_x[j][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*S[i][j];
      }
      
      /*** prior mean ***/
      mu_x[0] = 0;
      for (m = 0; m < *n_dim; m++)
	mu_x[0] += Z[i][m] * delta[m];
      
      bNormalReg2(D_x, x_i, var_epsilon, *n_pol, 1, 1, 1, mu_x, A0_x,
		 1, 1, 1, 1);

      /*** packing sampled x ***/
      X[i] = x_i[0];

      R_FlushConsole();
      R_CheckUserInterrupt();
    }

    /*** storing sampled x ***/
    if(main_loop > *burn) {
      if(keep == *thin) {
	if (*x_sd){
	  varx = var(X, *n_samp, 1);
	  xStore[ix++] = sqrt(varx);
	} else {
	  for (i = 0; i < *n_samp; i++)
	    xStore[ix++] = X[i];
	}
      }
    }
    /** end of sampling x **/

   /** Start sampling delta **/
    for (i = 0; i < *n_samp; i++) {
      for (m = 0; m < *n_dim; m++)
	D_delta[i][m] = Z[i][m];
      
      D_delta[i][*n_dim] = X[i];
    }

    bNormalReg2(D_delta, delta, sig2_x, *n_samp, *n_dim, 1, 1,
		mu_delta, A0_delta, 1, 1, 1, 1);

    if (main_loop > *burn) {
      if (keep == *thin) {
	for (m = 0; m < *n_dim; m++)
	  deltaStore[idelta++] = delta[m];
      }
    }
    /** end of sampling delta **/

    /** start sampling lambda and omega2 **/
    for (j = 0; j < *n_pol; j++) {
      for (k = 0; k < *n_act; k++) {

	/*** # of observations with treatment k for question j ***/
	n = 0;
	itemp = 0;
	for (i = 0; i < *n_samp; i++) {
	  if (k+1 == T[i][j]) {
	    stemp[n] = S[i][j];

	    for (m = 0; m < *n_dim; m++)
	      ztemp[itemp++] = Z[i][m];

	    n++;
	  }
	}
	
	if (n == 0)
	  continue;
	
	double **D_lambda = doubleMatrix(n + *n_dim, *n_dim+1);

	itemp = 0;
	for (l = 0; l < n; l++) {

	  for (m = 0; m < *n_dim; m++)
	    D_lambda[l][m] = ztemp[itemp++];

	  D_lambda[l][*n_dim] = stemp[l];
	}

	for (m = 0; m < *n_dim; m++) {
	  mu_lambda[m] = Theta[k][m];
	  A0_lambda[m][m] = (1 / phi2[(k * *n_dim) + m]);
	}
	
	bNormalReg2(D_lambda, lambda_jk, omega2_jk, n, *n_dim, 1, 1, mu_lambda,
		   A0_lambda, 1, *s0_omega2, *nu0_omega2, 0);

	for (m = 0; m < *n_dim; m++)
	  Lambda[j][k][m] = lambda_jk[m];

	Omega2[j][k] = omega2_jk[0];

	FreeMatrix(D_lambda, n + *n_dim);

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    } 

    /*** storing sampled lambda ***/
    if (main_loop > *burn) {
      if (keep == *thin) {
	for (j = 0; j < *n_pol; j++)
	  for (k = 0; k < *n_act; k++)
	    for (m = 0; m < *n_dim; m++)
	      lambdaStore[ilambda++] = Lambda[j][k][m];
      }
    }
    /** end of sampling lambda and omega2 **/

    /*** LIST EXPERIMENT MODEL SAMPLING BEGINS HERE ***/
    for (j = 0; j < *n_act; j++)
      treatSum[j] = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar for treated units */
      if (treat[i] > 0) {
	if (Ylist[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0;  Xpsi[i] = 0;  
	  for (j = 0; j < *n_dim; j++) {
	    Xdelta[i] += Z[i][j] * Theta[treat[i]-1][j];
	    Xpsi[i] += Z[i][j] * psi[j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Ylist[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(pt(Xdelta[i], *df, 1, 1) + dbinom(*J-1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp2 = exp(pt(Xdelta[i], *df, 1, 1) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp3 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Ylist[i] == 1)) { /* floor effects */
	    dtemp1 = exp(pt(Xdelta[i], *df, 1, 1) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp3 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    dtemp1 = exp(pt(Xdelta[i], *df, 1, 1) + 
			 dbinom(Ylist[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(Ylist[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Ylist[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Ylist[i];
	    }
	  }
	}
	treatSum[treat[i]-1]++;
      }
    }

    /* Sample psi */
    BinomLogit(Y0, Z, psi, *n_samp, *J, *n_dim, psi0, A0psi, psiPro, 1, psiCounter);
    if (main_loop > *burn) {
      if (keep == *thin) {
	for (j = 0; j < *n_dim; j++)
	  allresults[itempS++] = psi[j];
	allresults[itempS++] =  ((double) *psiCounter / (double) main_loop);
      }
    }

    /*** JOINT SAMPLING OF THETA ***/
    /** start sampling theta and phi2 **/
    for (k = 0; k < *n_act; k++) {

      /* sampling latent variable for list experiment */
      for (i = 0; i < treatSum[k]; i++) {
	dtemp = 0; 
	for (m = 0; m < *n_dim; m++) 
	  dtemp += Xtemp[k][i][m]*Theta[k][m];
	if (Zstar[k][i] == 0)
	  Xtemp[k][i][*n_dim] = TruncT(dtemp-1000,0,dtemp,*df,1,1);
	else
	  Xtemp[k][i][*n_dim] = TruncT(0,dtemp+1000,dtemp,*df,1,1);
	Prec[k][i] = rgamma(((double)*df + 1)/2, 
			    ((double)*df + (Xtemp[k][i][*n_dim] - dtemp)*(Xtemp[k][i][*n_dim] - dtemp))/2);
      }

      itemp = treatSum[k];
      for (m = 0; m < *n_dim; m++) {
	dtemp = 0;
	for (j = 0; j < *n_pol; j++) {
	  Xtemp[k][itemp][*n_dim] = Lambda[j][k][m];
	  dtemp += (Lambda[j][k][m] - Theta[k][m])*(Lambda[j][k][m] - Theta[k][m]);
	  for (i = 0; i < *n_dim; i++) {
	    if (i == m) 
	      Xtemp[k][itemp][i] = 1;
	    else
	      Xtemp[k][itemp][i] = 0;
	  }
	  itemp++;
	}
	phi2[(k * *n_dim) + m] = 1/rgamma(((double)*df + *n_pol)/2, ((double)*df + dtemp)/2);
	for (j = 1; j <= *n_pol; j++) {
	  Prec[k][itemp-j] = 1/phi2[(k * *n_dim) + m];
	}
      }

      bNormalReg1(Xtemp[k], Theta[k], treatSum[k] + *n_pol * *n_dim, *n_dim, 
		  1, dmu_theta, A0_theta, Prec[k], *df, 1);

      /*** storing sampled theta ***/
      if (main_loop > *burn) {
	if (keep == *thin) {
	  for (m = 0; m < *n_dim; m++) {
	    thetaStore[itheta++] = Theta[k][m];
	  }
	}
	
      }
      R_FlushConsole();
      R_CheckUserInterrupt();
    } /** end of sampling theta and phi2 **/

 
    /** update thinning counter **/
    if (keep == *thin) {
      keep = 1;
    } else {
      keep++;
    }

    /* printing */
    if (main_loop % 100 == 0) {
      Rprintf("Acceptance rate for psi: %4g\n", (j + 1), fprec((double) *psiCounter / (double) main_loop, 3));
    }

    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  Rprintf("End Gibbs Sampler\n");

  PutRNGstate();


  /* freeing memory */
  free(Y_j);
  free(beta);
  free(mu_beta);
  free(s);
  free(var_epsilon);
  free(mu_s);
  free(x_i);
  free(mu_x);
  free(lambda_jk);
  free(omega2_jk);
  free(mu_lambda);
  free(stemp);
  free(ztemp);
  free(theta_km);
  free(mu_theta);
  free(sig2_x);
  free(accept);
  free(temp_accept);

  FreeintMatrix(T, *n_samp);
  FreeintMatrix(Y, *n_samp);
  FreeMatrix(Z, *n_samp);
  FreeMatrix(S, *n_samp);
  FreeMatrix(Beta, *n_pol);
  FreeMatrix(Tau, *n_pol);
  FreeMatrix(Theta, *n_act);
  Free3DMatrix(Lambda, *n_pol, *n_act);
  FreeMatrix(Omega2, *n_pol);
  FreeMatrix(Mu_beta, *n_pol);
  FreeMatrix(A0_beta, 2);
  FreeMatrix(Ystar, *n_samp);
  FreeMatrix(U, *n_samp+2);
  FreeMatrix(D_s, 2);
  FreeMatrix(A0_s, 1);
  FreeMatrix(D_x, *n_pol+1);
  FreeMatrix(A0_x, 1);
  FreeMatrix(A0_lambda, 1);
  FreeMatrix(D_theta, *n_pol+1);
  FreeMatrix(A0_theta, *n_dim);
  FreeMatrix(A0_delta, *n_dim);
  FreeMatrix(D_delta, *n_samp + *n_dim);

  /* list experiment stuff */
  free(treatSum);
  free(psiCounter);
  free(Y0);
  free(Xdelta);
  free(Xpsi);
  FreeMatrix(A0psi, *n_dim);
  FreeMatrix(psiPro, *n_dim);
  FreeintMatrix(Zstar, *n_act);
  Free3DMatrix(Xtemp, *n_act, *n_samp + *n_pol * *n_dim);
  FreeMatrix(Prec, *n_act);
  
}

/* 
   Bayesian Robit constrained ictreg 

 */

void ictregBinomMultiRobit(int *Y,             /* outcome vector */
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
			   double *psiVar,     /* proposal variance for psi */
			   int *df,             /* degrees of freedom for t-distribution */
			   int *ceiling,       /* ceiling effects */
			   int *floor,         /* floor effects */
			   int *burnin,        /* number of burnins */
			   int *keep,          /* keep every *th draw */
			   int *verbose,       /* want to print progress? */
			   double *allresults  /* storage for results */
			   ) {

  int i, j, k, main_loop, itemp, itempS, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int *psiCounter = intArray(1);        /* acceptance ratio for psi */
  int *treatSum = intArray(*tmax);
  double dtemp, dtemp1, dtemp2, dtemp3;

  /** get random seed **/
  GetRNGstate();

  /** Parameters for sensitive items **/
  int n_dim = *n_cov;

  double **deltaMatrix = doubleMatrix(*tmax, n_dim);

  itemp = 0; 
  for (i = 0; i < *tmax; i++) 
    for (j = 0; j < n_dim; j++)
      deltaMatrix[i][j] = delta[itemp++];

  /** Data **/
  int **Zstar = intMatrix(*tmax, *n_samp);
  int *Y0 = intArray(*n_samp);
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim); 
  double ***Xtemp = doubleMatrix3D(*tmax, *n_samp, n_dim); 
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  for (i = 0; i < *tmax; i++)
    treatSum[i] = 0;

  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Y[i];
    if (treat[i] > 0) {
      for (j = 0; j < *n_cov; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = X[i][j];
      treatSum[treat[i]-1]++;
      if ((ceiling[treat[i]-1] == 1) && (Y[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Y[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
  }

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

  /** Proposal precisoin **/
  double **psiPro = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      psiPro[i][j] = psiVar[itemp++];

  /** MCMC **/
  itempS = 0; psiCounter[0] = 0; 
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... *** \n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *tmax; j++)
      treatSum[j] = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar for treated units */
      if (treat[i] > 0) {
	if (Y[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0;  Xpsi[i] = 0;  
	  for (j = 0; j < *n_cov; j++) {
	    Xdelta[i] += X[i][j] * deltaMatrix[treat[i]-1][j];
	    Xpsi[i] += X[i][j] * psi[j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Y[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(pt(Xdelta[i], *df, 1, 1) + dbinom(*J-1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp2 = exp(pt(Xdelta[i], *df, 1, 1) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp3 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Y[i] == 1)) { /* floor effects */
	    dtemp1 = exp(pt(Xdelta[i], *df, 1, 1) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp3 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    dtemp1 = exp(pt(Xdelta[i], *df, 1, 1) + 
			 dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(pt(Xdelta[i], *df, 0, 1) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Y[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Y[i];
	    }
	  }
	}
	treatSum[treat[i]-1]++;
      }
    }

    /* Sample delta */
    for (k = 0; k < *tmax; k++) {
      RobitGibbs(Zstar[k], Xtemp[k], deltaMatrix[k], treatSum[k],
		 n_dim, 1, m0delta[k], A0delta[k], *df, 1);
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
	allresults[itempS++] = ((double) *psiCounter / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("\n    Control items model: %3g\n", 
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
  free(treatSum);
  free(Y0);
  free(Xdelta);
  free(Xpsi);
  FreeMatrix(m0delta, *tmax);
  FreeMatrix(deltaMatrix, *tmax);
  FreeintMatrix(Zstar, *tmax);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(Xtemp, *tmax, *n_samp);
  Free3DMatrix(A0delta, *tmax, n_dim);
  FreeMatrix(A0psi, *n_cov);
  FreeMatrix(psiPro, *n_cov);
}


