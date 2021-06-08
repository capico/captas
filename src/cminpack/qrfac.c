#include "cminpack.h"


void qrfac(int m, int n, double **a, int pivot, int *ipvt, double *rdiag,
		   double *acnorm, double *wa)
/**
    subroutine qrfac

    this subroutine uses householder transformations with column
    pivoting (optional) to compute a qr factorization of the
    m by n matrix a. that is, qrfac determines an orthogonal
    matrix q, a permutation matrix p, and an upper trapezoidal
    matrix r with diagonal elements of nonincreasing magnitude,
    such that a*p = q*r. the householder transformation for
    column k, k = 1,2,...,min(m,n), is of the form

                       t
       i - (1/u(k))*u*u

    where u has zeros in the first k-1 positions. the form of
    this transformation and the method of pivoting first
    appeared in the corresponding linpack subroutine.

    the subroutine statement is pivot,ipvt,lipvt,rdiag,acnorm,wa)

    where

    m is a positive integer input variable set to the number
     of rows of a.

    n is a positive integer input variable set to the number
     of columns of a.

    a is an m by n array. on input a contains the matrix for
     which the qr factorization is to be computed. on output
     the strict upper trapezoidal part of a contains the strict
     upper trapezoidal part of r, and the lower trapezoidal
     part of a contains a factored form of q (the non-trivial
     elements of the u vectors described above).

    lda is a positive integer input variable not less than m
     which specifies the leading dimension of the array a.

    pivot is a logical input variable. if pivot is set true,
     then column pivoting is enforced. if pivot is set false,
     then no column pivoting is done.

    ipvt is an integer output array of length lipvt. ipvt
     defines the permutation matrix p such that a*p = q*r.
     column j of p is column ipvt(j) of the identity matrix.
     if pivot is false, ipvt is not referenced.

    lipvt is a positive integer input variable. if pivot is false,
     then lipvt may be as small as 1. if pivot is true, then
     lipvt must be at least n.

    rdiag is an output array of length n which contains the
     diagonal elements of r.

    acnorm is an output array of length n which contains the
     norms of the corresponding columns of the input matrix a.
     if this information is not needed, then acnorm can coincide
     with rdiag.

    wa is a work array of length n. if pivot is false, then wa
     can coincide with rdiag.

    subprograms called

    minpack-supplied ... dpmpar,enorm

    fortran-supplied ... dmax1,dsqrt,min0

    argonne national laboratory. minpack project. march 1980.
    burton s. garbow, kenneth e. hillstrom, jorge j. more
*/
{
    int i, j, jp1, k, kmax, minmn;
    double ajnorm, epsmch, sum, temp;

	/* get machine precision */
    epsmch = MACHEPS;
	/* compute the initial column norms and initialize several arrays */
    for(j = 0; j < n; j++){
        acnorm[j] = enorm(m, &a[j][0]);
        rdiag[j]  = acnorm[j];
        wa[j]     = rdiag[j];
        if(pivot)
			ipvt[j] = j;
	}
	/* reduce a to r with householder transformations */
    minmn = (m < n) ? m : n;
    for(j = 0; j < minmn; j++){
        if(pivot){
			/* bring column with largest norm into the pivot position */
            kmax = j;
            for(k = j; k < n; k++) /* +++++++++++++++++ */
                if(rdiag[k] > rdiag[kmax])
					kmax = k;
            if(kmax != j){
                for(i = 0; i < m; i++){
                    temp       = a[j][i];
                    a[j][i]    = a[kmax][i];
                    a[kmax][i] = temp;
				}
                rdiag[kmax] = rdiag[j];
                wa[kmax]    = wa[j];
                k           = ipvt[j];
                ipvt[j]     = ipvt[kmax];
                ipvt[kmax]  = k;
            }
        }
		/* compute the householder transformation to reduce the
		j-th column of a to a multiple of the j-th unit vector. */
		ajnorm = enorm(m - j, &a[j][j]);
		if(ajnorm != 0.0){
            if(a[j][j] < 0.0)
				ajnorm = -ajnorm;
            for(i = j; i < m; i++)
                a[j][i] /= ajnorm;
            a[j][j] += 1.0;
			/* apply the transformation to the remaining columns and
			update the norms. */
            jp1 = j + 1;
            if(n > jp1){
                for(k = jp1; k < n; k++){
                    sum = 0.0;
                    for(i = j; i < m; i++)
                        sum += a[j][i]*a[k][i];
                    temp = sum / a[j][j];
                    for(i = j; i < m; i++)
                        a[k][i] -=temp*a[j][i];
                    if(!pivot || !rdiag[k])
						continue;
                    temp = a[k][j] / rdiag[k];
                    rdiag[k] *= sqrt(max(0.0, 1.0 - temp*temp));
                    if(0.5 * (rdiag[k] * rdiag[k] / (wa[k] * wa[k])) > epsmch)
						continue;
					rdiag[k] = enorm(m - j - 1, &a[k][jp1]);
                    wa[k]    = rdiag[k];
                }
            }
        }
        rdiag[j] = -ajnorm;
    }

	return;
}
