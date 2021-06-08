#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <malloc.h>

#ifndef CMINPACK_H
#define CMINPACK_H

#define MACHEPS DBL_EPSILON
#define MINMAG  DBL_MIN
#define MAXMAG  DBL_MAX

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)	(((a) < (b)) ? (a) : (b))


/* This structure must contain at least two vectors with the abscissas and
ordinates of the data to be adjusted. Other vectors may contain the weights of
the observations or other values that are convenient when used by other
functions. */
typedef struct {
    double *t;  /* independent variable */
    double *p;  /* observations */

	/* weights or temporary arrays */
	double *dp;
	double *dt;
	double *dte;
	double *derl;
	double *derg;
} lmdatatype;   /* structure for the data to be fitted */


double enorm(int n, double *x);

void lmpar(int n, double **fjac, int *ipvt, double *diag, double *qtf,
           double delta, double *par, double *wa1, double *wa2, double *wa3,
           double *wa4);

void qrfac(int m, int n, double **a, int pivot, int *ipvt, double *rdiag,
           double *acnorm, double *wa);

void qrsolv(int n, double **r, int *ipvt, double *diag, double *qtb, double *x,
			double *sdiag, double *wa);

void lmder(void f(), lmdatatype *data, void *modelpar, int m, int n, double *x,
		   double *fvec, double ftol, double xtol, double gtol, int maxfev,
           double *diag, int mode, double factor, int *info, int *nfev, int *njev,
           double **fjac, int *ipvt, double *qtf, double *wa1, double *wa2,
           double *wa3, double *wa4);

void lmder1(void fcn(), lmdatatype *data, void *modelpar, int m, int n, double *x,
		   double *fvec, double tol, int *info, int *nfev, int *njev);

void covar(int n, double **r, int *ipvt, double tol, double *wa);

#endif
