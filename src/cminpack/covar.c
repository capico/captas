#include "cminpack.h"


void covar(int n, double **r, int *ipvt, double tol, double *wa)
/**
    subroutine covar

    given an m by n matrix a, the problem is to determine
    the covariance matrix corresponding to a, defined as

                   t
          inverse(a *a) .

    this subroutine completes the solution of the problem
    if it is provided with the necessary information from the
    qr factorization, with column pivoting, of a. that is, if
    a*p = q*r, where p is a permutation matrix, q has orthogonal
    columns, and r is an upper triangular matrix with diagonal
    elements of nonincreasing magnitude, then covar expects
    the full upper triangle of r and the permutation matrix p.
    the covariance matrix is then computed as

                     t     t
          p*inverse(r *r)*p  .

    if a is nearly rank deficient, it may be desirable to compute
    the covariance matrix corresponding to the linearly independent
    columns of a. to define the numerical rank of a, covar uses
    the tolerance tol. if l is the largest integer such that

          abs(r(l,l)) .gt. tol*abs(r(1,1)) ,

    then covar computes the covariance matrix corresponding to
    the first l columns of r. for k greater than l, column
    and row ipvt(k) of the covariance matrix are set to zero.

    the subroutine statement is

      subroutine covar(n,r,ldr,ipvt,tol,wa)

    where

      n is a positive integer input variable set to the order of r.

      r is an n by n array. on input the full upper triangle must
        contain the full upper triangle of the matrix r. on output
        r contains the square symmetric covariance matrix.

      ldr is a positive integer input variable not less than n
        which specifies the leading dimension of the array r.

      ipvt is an integer input array of length n which defines the
        permutation matrix p such that a*p = q*r. column j of p
        is column ipvt(j) of the identity matrix.

      tol is a nonnegative input variable used to define the
        numerical rank of a in the manner described above.

      wa is a work array of length n.

    subprograms called

      fortran-supplied ... dabs

    argonne national laboratory. minpack project. august 1980.
    burton s. garbow, kenneth e. hillstrom, jorge j. more
*/
{
	double tolr, temp;
	int i, ii, j, jj, k, km1, l = -1, sign;

	/* form the inverse of r in the full upper triangle of r. */
	tolr = tol*fabs(r[0][0]);
	for(k = 0; k < n; k++){

		if(fabs(r[k][k]) <= tolr)
			break;

		r[k][k] = 1.0/r[k][k];
		km1 = k - 1;

		if(km1 < 0){

			l = k;
		}
		else{

			for(j = 0; j <= km1; j++){
				temp = r[k][k]*r[k][j];
				r[k][j] = 0.0;

				for(i = 0; i <= j; i++){
					r[k][i] = r[k][i] - temp*r[j][i];
				}
			}

			l = k;
		}
	}

	/*************************************************************************/
	/* form the full upper triangle of the inverse of (r transpose)*r in the
	full upper triangle of r. */
	if(l >= 0){

		for(k = 0; k <= l; k++){
			km1 = k - 1;

			if(km1 >= 0){

				for(j = 0; j <= km1; j++){
					temp = r[k][j];

					for(i = 0; i <= j; i++){
						r[j][i] = r[j][i] + temp*r[k][i];
					}
				}
			}

			temp = r[k][k];

			for(i = 0; i <= k; i++){
				r[k][i] = temp*r[k][i];
			}
		}
	}

	/*************************************************************************/
	/* form the full lower triangle of the covariance matrix in the strict
	lower triangle of r and in wa. */
	for(j = 0; j < n; j++){

		jj = ipvt[j];
		sign = (j > l);

		for(i = 0; i < j; i++){

			if(sign)
				r[j][i] = 0.0;

			ii = ipvt[i];

			if(ii > jj)
				r[jj][ii] = r[j][i];

			if(ii < jj)
				r[ii][jj] = r[j][i];
		}

		wa[jj] = r[j][j];
	}

	/*************************************************************************/
	/* symmetrize the covariance matrix in r. */
	for(j = 0; j < n; j++){

		for(i = 0; i < j; i++){
			r[j][i] = r[i][j];
		}

		r[j][j] = wa[j];
	}

	return;
}
