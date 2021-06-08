#include "cminpack.h"


void lmpar(int n, double **r, int ipvt[], double diag[], double qtb[],
		   double delta, double *par, double x[], double sdiag[], double wa1[],
		   double wa2[])
/**
    subroutine lmpar

    given an m by n matrix a, an n by n nonsingular diagonal
    matrix d, an m-vector b, and a positive number delta,
    the problem is to determine a value for the parameter
    par such that if x solves the system

    a*x = b,     sqrt(par)
    *d*x = 0,

    in the least squares sense, and dxnorm is the euclidean
    norm of d*x, then either par is zero and

    (dxnorm-delta) .le. 0.1*delta,

    or par is positive and

    abs(dxnorm-delta) .le. 0.1*delta .

    this subroutine completes the solution of the problem
    if it is provided with the necessary information from the
    qr factorization, with column pivoting, of a. that is, if
    a*p = q*r, where p is a permutation matrix, q has orthogonal
    columns, and r is an upper triangular matrix with diagonal
    elements of nonincreasing magnitude, then lmpar expects
    the full upper triangle of r, the permutation matrix p,
    and the first n components of (q transpose)
    *b. on output
    lmpar also provides an upper triangular matrix s such that

    t   t                   t
    p *(a *a + par*d*d)*p = s *s .

    s is employed within lmpar and may be of separate interest.

    only a few iterations are generally needed for convergence
    of the algorithm. if, however, the limit of 10 iterations
    is reached, then the output par will contain the best
    value obtained so far.

    the subroutine statement is

    subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
    wa1,wa2)

    where

    n is a positive integer input variable set to the order of r.

    r is an n by n array. on input the full upper triangle
    must contain the full upper triangle of the matrix r.
    on output the full upper triangle is unaltered, and the
    strict lower triangle contains the strict upper triangle
    (transposed) of the upper triangular matrix s.

    ldr is a positive integer input variable not less than n
    which specifies the leading dimension of the array r.

    ipvt is an integer input array of length n which defines the
    permutation matrix p such that a*p = q*r. column j of p
    is column ipvt(j) of the identity matrix.

    diag is an input array of length n which must contain the
    diagonal elements of the matrix d.

    qtb is an input array of length n which must contain the first
    n elements of the vector (q transpose)*b.

    delta is a positive input variable which specifies an upper
    bound on the euclidean norm of d*x.

    par is a nonnegative variable. on input par contains an
    initial estimate of the levenberg-marquardt parameter.
    on output par contains the final estimate.

    x is an output array of length n which contains the least
    squares solution of the system a*x = b, sqrt(par)*d*x = 0,
    for the output par.

    sdiag is an output array of length n which contains the
    diagonal elements of the upper triangular matrix s.

    wa1 and wa2 are work arrays of length n.

    subprograms called

    minpack-supplied ... dpmpar,enorm,qrsolv

    fortran-supplied ... dabs,dmax1,dmin1,dsqrt

    argonne national laboratory. minpack project. march 1980.
    burton s. garbow, kenneth e. hillstrom, jorge j. more
*/
{
    int i, iter, j, jm1, jp1, k, l, nsing;
    double dxnorm, dwarf, fp, gnorm, parc, parl, paru;
    double sum, temp;

    dwarf = MINMAG;
	/* compute and store in x the gauss-newton direction. if the
	jacobian is rank-deficient, obtain a least squares solution. */
    nsing = n;
    for(j = 0; j < n; j++){
        wa1[j] = qtb[j];

		if((r[j][j] == 0.0) && (nsing == n))
            nsing = j;

		if(nsing < n) wa1[j] = 0.0;
    }
    if(nsing >= 1){
        for(k = 0; k < nsing; k++){
            j = nsing - k - 1;
            wa1[j] /= r[j][j];
            temp    = wa1[j];

			if(j < 1)
				continue;

			for(i = 0; i < j; i++)
                wa1[i] -= r[j][i] * temp;
        }
    }
    for(j = 0; j < n; j++){
        l    = ipvt[j];
        x[l] = wa1[j];
    }
    /* initialize the iteration counter.
    evaluate the function at the origin, and test
    for acceptance of the gauss-newton direction. */
    iter = 0;
    for(j = 0; j < n; j++){
        wa2[j] = diag[j] * x[j];
    }

	dxnorm = enorm(n, wa2);
    fp     = dxnorm - delta;

	if(fp <= 0.1*delta){
        if(iter == 0)
            *par = 0.0;
        return;
    }
	/* if the jacobian is not rank deficient, the newton
	step provides a lower bound, parl, for the zero of
	the function. otherwise set this bound to zero. */
    parl = 0.0;
    if(nsing >= n){
        for(j = 0; j < n; j++){
            l = ipvt[j];
            wa1[j] = diag[l] * wa2[l] / dxnorm;
        }

		for(j = 0; j < n; j++){
            sum = 0.0;
			if(j >= 1){
				for(i = 0; i < j; i++)
					sum += r[j][i] * wa1[i];
            }
            wa1[j] = (wa1[j] - sum) / r[j][j];
        }
        temp = enorm(n,wa1);
        parl = ((fp / delta) / temp) / temp;
    }

	/* calculate an upper bound, paru, for the zero of the function. */
    for(j = 0;j < n; j++){
		sum = 0.0;
		for(i = 0; i <= j; i++)
			sum += r[j][i] * qtb[i];

		l      = ipvt[j];
		wa1[j] = sum / diag[l];
    }
    gnorm = enorm(n,wa1);
    paru  = gnorm / delta;

	if(paru == 0.0)
		paru = dwarf / min(delta,0.1);
	/* if the input par lies outside of the interval (parl,paru),
	set par to the closer endpoint. */
    *par = max(*par, parl);
    *par = min(*par, paru);
    if(*par == 0.0)
		*par = gnorm / dxnorm;
	/* beginning of an iteration. */
    while(1){
        iter++;
		/* evaluate the function at the current value of par. */
        if(*par == 0.0)
			*par = max(dwarf, 0.001 * paru);

		temp = sqrt(*par);
        for(j = 0; j < n; j++)
            wa1[j] = temp * diag[j];

		qrsolv(n, r, ipvt, wa1, qtb, x, sdiag, wa2);
        for(j = 0; j < n; j++)
            wa2[j] = diag[j] * x[j];

		dxnorm = enorm(n, wa2);
        temp   = fp;
        fp     = dxnorm - delta;
		/*  if the function is small enough, accept the current value
		of par. also test for the exceptional cases where parl
		is zero or the number of iterations has reached 10. */
        if( (fabs(fp) <= 0.1*delta) || (parl == 0.0) && (fp <= temp)
            && (temp > 0.0) || (iter == 10) ){
				if(iter == 0)
					*par = 0.0;
				return;
        }

		/* compute the newton correction. */
        for(j = 0; j < n; j++){
            l      = ipvt[j];
            wa1[j] = diag[l] * wa2[l] / dxnorm;
        }
        for(j = 0; j < n; j++){
            wa1[j] /= sdiag[j];
            temp    = wa1[j];
            jp1     = j + 1;

			if(jp1 < n)
				for(i = jp1; i < n; i++)
					wa1[i] -= r[j][i] * temp;
        }
        temp = enorm(n, wa1);
        parc = ((fp/delta) / temp) / temp;
		/* depending on the sign of the function, update parl or paru. */
        if(fp > 0.0)
            parl = max(parl, *par);
        if(fp < 0.0)
            paru = min(paru, *par);
		/* compute an improved estimate for par. */
        *par = max(parl, *par+parc);
    }

	return;
}
