#include "cminpack.h"


void lmder(void f(),
		   lmdatatype *data, void *modelpar, int m, int n, double *x,
		   double *fvec, double ftol, double xtol, double gtol, int maxfev,
		   double *diag, int mode, double factor, int *info, int *nfev, int *njev,
		   double **fjac, int *ipvt, double *qtf, double *wa1, double *wa2,
		   double *wa3, double *wa4)

/**
   subroutine lmder

   the purpose of lmder is to minimize the sum of the squares of
   m nonlinear functions in n variables by a modification of
   the levenberg-marquardt algorithm. the user must provide a
   subroutine which calculates the functions and the jacobian.

   the subroutine statement is

   subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
   maxfev,diag,mode,factor,nprint,info,nfev,
   njev,ipvt,qtf,wa1,wa2,wa3,wa4)

   where

   fcn is the name of the user-supplied subroutine which
   calculates the functions and the jacobian. fcn must
   be declared in an external statement in the user
   calling program, and should be written as follows.

   subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
   integer m,n,ldfjac,iflag
   double precision x(n),fvec(m),fjac(ldfjac,n)
   ----------
   if iflag = 1 calculate the functions at x and
   return this vector in fvec. do not alter fjac.
   if iflag = 2 calculate the jacobian at x and
   return this matrix in fjac. do not alter fvec.
   ----------
   return
   end

   the value of iflag should not be changed by fcn unless
   the user wants to terminate execution of lmder.
   in this case set iflag to a negative integer.

   m is a positive integer input variable set to the number
   of functions.

   n is a positive integer input variable set to the number
   of variables. n must not exceed m.

   x is an array of length n. on input x must contain
   an initial estimate of the solution vector. on output x
   contains the final estimate of the solution vector.

   fvec is an output array of length m which contains
   the functions evaluated at the output x.

   fjac is an output m by n array. the upper n by n submatrix
   of fjac contains an upper triangular matrix r with
   diagonal elements of nonincreasing magnitude such that

   t     t           t
   p *(jac *jac)*p = r *r,

   where p is a permutation matrix and jac is the final
   calculated jacobian. column j of p is column ipvt(j)
   (see below) of the identity matrix. the lower trapezoidal
   part of fjac contains information generated during
   the computation of r.

   ldfjac is a positive integer input variable not less than m
   which specifies the leading dimension of the array fjac.

   ftol is a nonnegative input variable. termination
   occurs when both the actual and predicted relative
   reductions in the sum of squares are at most ftol.
   therefore, ftol measures the relative error desired
   in the sum of squares.

   xtol is a nonnegative input variable. termination
   occurs when the relative error between two consecutive
   iterates is at most xtol. therefore, xtol measures the
   relative error desired in the approximate solution.

   gtol is a nonnegative input variable. termination
   occurs when the cosine of the angle between fvec and
   any column of the jacobian is at most gtol in absolute
   value. therefore, gtol measures the orthogonality
   desired between the function vector and the columns
   of the jacobian.

   maxfev is a positive integer input variable. termination
   occurs when the number of calls to fcn with iflag = 1
   has reached maxfev.

   diag is an array of length n. if mode = 1 (see
   below), diag is internally set. if mode = 2, diag
   must contain positive entries that serve as
   multiplicative scale factors for the variables.

   mode is an integer input variable. if mode = 1, the
   variables will be scaled internally. if mode = 2,
   the scaling is specified by the input diag. other
   values of mode are equivalent to mode = 1.

   factor is a positive input variable used in determining the
   initial step bound. this bound is set to the product of
   factor and the euclidean norm of diag*x if nonzero, or else
   to factor itself. in most cases factor should lie in the
   interval (.1,100.).100. is a generally recommended value.

   nprint is an integer input variable that enables controlled
   printing of iterates if it is positive. in this case,
   fcn is called with iflag = 0 at the beginning of the first
   iteration and every nprint iterations thereafter and
   immediately prior to return, with x, fvec, and fjac
   available for printing. fvec and fjac should not be
   altered. if nprint is not positive, no special calls
   of fcn with iflag = 0 are made.

   info is an integer output variable. if the user has
   terminated execution, info is set to the (negative)
   value of iflag. see description of fcn. otherwise,
   info is set as follows.

   info = 0  improper input parameters.

   info = 1  both actual and predicted relative reductions
   in the sum of squares are at most ftol.

   info = 2  relative error between two consecutive iterates
   is at most xtol.

   info = 3  conditions for info = 1 and info = 2 both hold.

   info = 4  the cosine of the angle between fvec and any
   column of the jacobian is at most gtol in
   absolute value.

   info = 5  number of calls to fcn with iflag = 1 has
   reached maxfev.

   info = 6  ftol is too small. no further reduction in
   the sum of squares is possible.

   info = 7  xtol is too small. no further improvement in
   the approximate solution x is possible.

   info = 8  gtol is too small. fvec is orthogonal to the
   columns of the jacobian to machine precision.

   nfev is an integer output variable set to the number of
   calls to fcn with iflag = 1.

   njev is an integer output variable set to the number of
   calls to fcn with iflag = 2.

   ipvt is an integer output array of length n. ipvt
   defines a permutation matrix p such that jac*p = q*r,
   where jac is the final calculated jacobian, q is
   orthogonal (not stored), and r is upper triangular
   with diagonal elements of nonincreasing magnitude.
   column j of p is column ipvt(j) of the identity matrix.

   qtf is an output array of length n which contains
   the first n elements of the vector (q transpose)*fvec.

   wa1, wa2, and wa3 are work arrays of length n.

   wa4 is a work array of length m.

   subprograms called

   user-supplied ...... fcn

   minpack-supplied ... dpmpar,enorm,lmpar,qrfac

   fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

   argonne national laboratory. minpack project. march 1980.
   burton s. garbow, kenneth e. hillstrom, jorge j. more
*/
{
    int i, iflag, iter, j, l;
    double actred, delta = 1.0, dirder, epsmch, fnorm, fnorm1, gnorm;
    double par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm = 0.0;

	/* initialize */
	/* epsmch is the machine precision. */
    epsmch = MACHEPS;

    *info = 0;
    iflag = 0;
    *nfev = 0;
	*njev = 0;

	/* check for input parameter errors */
    if((n <= 0) || (m < n) || (maxfev <= 0) || (ftol < 0.0) || (xtol < 0.0)
		|| (gtol < 0.0) || (maxfev <= 0) || (factor <= 0.0))
		return;

	if(mode == 2) {
        for (j = 0; j < n; j++)
            if (diag[j] <= 0)
				return;
    }

	/* evaluate the function at the starting point and calculate its norm */
    iflag = 1;
    f(m, n, x, fvec, fjac, data, modelpar, &iflag);
    *nfev = 1;
    if(iflag < 0) {
        *info = iflag;
        return;
    }
    fnorm = enorm(m, fvec);

	/* initialize levenberg-marquardt parameter and iteration counter. */
    par  = 0;
    iter = 1;

	/* beginning of the outer loop */
    while(1){
		/* calculate jacobian matrix */
        iflag = 2;
        f(m, n, x, fvec, fjac, data, modelpar, &iflag);
        (*njev)++;
        if(iflag < 0) {
            *info = iflag;
            return;
        }

		/* if requested, call fcn to enable printing of iterates. */
		iflag = 0;
		/*printf("%f ", delta);*/
		//f(m, n, x, fvec, fjac, data, modelpar, &iflag); <- uncomment to print iterations
		/*printf("fnorm %.15e\n", enorm(m,fvec));
		for(i = 0; i < n; i++){
			printf("%.15e ", x[i]);
		}
		printf("\n");*/
		/*iflag = 0;
		if((iter%nprint) == 0)
			f(m, n, x, fvec, &iflag);
		if (iflag < 0) {
            *info = iflag;
            return;
        }*/


		/* compute the qr factorization of the jacobian */
        qrfac(m, n, fjac, 1, ipvt, wa1, wa2, wa3);

		/* on the first iteration and if mode is 1, scale according
		to the norms of the columns of the initial jacobian. */
        if(iter == 1) {
            if(mode != 2) {
                for(j = 0; j < n; j++) {
                    diag[j] = wa2[j];
                    if(wa2[j] == 0.0)
						diag[j] = 1.0;
                }
            }
			/*  on the first iteration, calculate the norm of the scaled x
			and initialize the step bound delta. */
            for(j = 0; j < n; j++)
                wa3[j] = diag[j] * x[j];
            xnorm = enorm(n, wa3);
            delta = factor * xnorm;
            if(delta == 0)
				delta = factor;
        }

		/* form (q transpose)*fvec and store the first n components in qtf. */
        for(i = 0; i < m; i++)
            wa4[i] = fvec[i];
        for(j = 0; j < n; j++){
            if(fjac[j][j] != 0.0){
                sum = 0.0;
                for(i = j; i < m; i++)
                    sum += fjac[j][i] * wa4[i];
                temp = -sum / fjac[j][j];
                for(i = j; i < m; i++)
                    wa4[i] += fjac[j][i] * temp;
            }
            fjac[j][j] = wa1[j];
            qtf[j]     = wa4[j];
        }

		/* compute the norm of the scaled gradient */
        gnorm = 0.0;
        if(fnorm != 0.0){
            for(j = 0; j < n; j++){
                l = ipvt[j];
                if(wa2[l] == 0.0)
					continue;
                sum = 0.0;
				/* ++++++++++++++++++++++++++++ */
                for(i = 0; i <= j; i++)
                    sum += fjac[j][i] * qtf[i] / fnorm;
                gnorm = max(gnorm, fabs(sum/wa2[l]));
            }
        }

		/* test for convergence of the gradient norm */
        if(gnorm <= gtol)
			*info = 4;
        if(*info != 0){
			/* ++++++++++++++ */
            if(*info < 0)
					*info = iflag;
			return;
        }

		/* rescale if necessary */
        if(mode != 2){
			for(j = 0; j < n; j++)
				diag[j] = max(diag[j], wa2[j]);
        }

		/* beginning of the inner loop */
        do{
			/* determine the levenberg-marquardt parameter */
            lmpar(n, fjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4);

			/* store the direction p and x + p. calculate the norm of p. */
			for(j = 0; j < n; j++){
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            pnorm = enorm(n, wa3);

			/* on the first iteration, adjust the initial step bound. */
            if(iter == 1)
				delta = 1.0*min(delta, pnorm); /* +++++++++++++++ */

			/* evaluate the function at x + p and calculate its norm. */
            iflag = 1;
            f(m, n, wa2, wa4, fjac, data, modelpar, &iflag);
            (*nfev)++;
            if(iflag < 0){
                *info = iflag;
                return;
            }
            fnorm1 = enorm(m, wa4);

			/* compute the scaled actual reduction. */
            actred = -1.0;
            if(0.1 * fnorm1 < fnorm)
                actred = 1.0 - (fnorm1*fnorm1/(fnorm*fnorm));

			/* compute the scaled predicted reduction and
			the scaled directional derivative. */
			for(j = 0; j < n; j++){
                wa3[j] = 0.0;
                l      = ipvt[j];
                temp   = wa1[l];
				/* ++++++++++++++++++++++ */
                for(i = 0; i <= j; i++)
                    wa3[i] += fjac[j][i] * temp;
            }
            temp1  = enorm(n, wa3) / fnorm;
            temp2  = sqrt(par) * pnorm / fnorm;
            prered = temp1*temp1 + 2.0*temp2*temp2;
            dirder = -(temp1*temp1 + temp2*temp2);

			/* compute the ratio of the actual to the predicted reduction. */
            ratio = 0.0;
            if(prered != 0.0)
				ratio = actred/prered;

			/* update the step bound. */
            if(ratio <= 0.25){
				/* ++++++++++++++++++++++ */
                if(actred >= 0.0)
					temp = 0.5;
                if(actred < 0.0)
					temp = 0.5*dirder/(dirder + 0.5*actred);
				/* +++++++++++++++++ */
				if((0.1*fnorm1 >= fnorm) || (temp < 0.1))
					temp = 0.1;
                delta = temp * min(delta, pnorm/0.1);
                par  /= temp;
            }
            else{
                if((par == 0.0) || (ratio >= 0.75)){
                    delta = pnorm / 0.5;
                    par  *= 0.5;
                }
            }

			/* test for successful iteration. */
            if(ratio >= 0.0001){
				/* successful iteration. update x, fvec, and their norms. */
                for(j = 0; j < n; j++){
					x[j]   = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }

				for(i = 0; i < m; i++)
                    fvec[i] = wa4[i];
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }

			/* tests for convergence. */
			*info = 0;
			if(fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1)
				*info = 1;
			if(delta <= xtol * xnorm)
				*info += 2;
			if(*info != 0)
				return;

			/* tests for termination and stringent tolerances. */
            if(*nfev >= maxfev)
				*info = 5;
            if((fabs(actred) <= epsmch) && (prered <= epsmch) &&
                (0.5*ratio <= 1.0))
				*info = 6;
            if(delta <= epsmch*xnorm)
				*info = 7;
            if(gnorm <= epsmch)
				*info = 8;
			/* +++++++++++++++++++++ */
            if(*info != 0) {
                if(*info < 0)
					*info = iflag;
                return;
            }
			/* end of the inner loop. repeat if iteration unsuccessful. */
        } while(ratio <= 0.0001);

		/* end of the outer loop. */
    }

	return;
}
