#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "stehfest.h"

/*---------------------------------------------------------------------------*/
/**
@brief    calculates Stehfest algorithm coefficients
@param    nstehfest	number of Stehfest terms
@param    ptr to current Stehfest coefficients vector (allocated memory)
or to NULL
@return   ptr to new Stehfest coefficients vector

This function returns a pointer to a dynamically allocated vector
containing the "v" coefficients for the Stehfest algorithm. If pointer *v is
not NULL, the function frees the vector and returns a new vector, with
nstehfest elements.

References:
Stehfest (1970a), Algorithm 368: Numerical inversion of Laplace transforms
[D5], Communications of the ACM, Volume 13 Issue 1, Jan. 1970, Pages 47-49

Stehfest (1970b), Remark on algorithm 368 [D5]: Numerical inversion of
Laplace transforms, Communications of the ACM, Volume 13 Issue 10, Oct. 1970,
Page 624
*/
/*---------------------------------------------------------------------------*/
double* stehfest_init(const int nstehfest, double *v)
{
	int i, k, n, nh, nv, kmax, sign = 1;
	double *g, *h, *vtmp;

	/* clear previous values */
	if(v != NULL){
		free(v);
	}

	/* Stehfest parameter */
	n = nstehfest;

	nh	= n / 2;
	nv	= 2 * nh;

	g = (double*)calloc(nv + 1, sizeof(double));
	h = (double*)calloc(nh + 1, sizeof(double));

	vtmp = (double*)calloc(nv, sizeof(double));

	g[0] = 1.0;
	for(i = 1; i <= nv; i++){
		g[i] = g[i-1]*i;
	}

	h[1] = 2/g[nh-1];
	for (i = 2; i <= nh; i++){
		h[i] = pow((double)i, nh)*g[2*i]/(g[nh - i]*g[i]*g[i - 1]);
	}

	if ((nh % 2) == 0)
		sign = -1.0;

	for(i = 1; i <= nv; i++){
		vtmp[i-1] = 0.0;
		kmax = (i < nh)? i: nh;
		for(k = (i + 1)/2; k <= kmax; k++)
			vtmp[i-1] += h[k] / (g[i-k]*g[2*k-i]);
		vtmp[i-1] *= sign;
		sign = -sign;
	}

	free(g);
	free(h);

	return vtmp;
}


/*---------------------------------------------------------------------------*/
/**
@brief    calculates Stehfest numerical inverse of Laplace Transforms
@param    f, pointer to the Laplace Transform function, f(par, u)
@param    par, pointer to the parameters of function f
@param    nstehfest	number of Stehfest terms
@param    ptr to current Stehfest coefficients vector
@param    t, time
@return   numerical inverse of f(par, u) at time t: F(par, t)

This function returns the numerical inverse of the Laplace Transform function
f(par, u). The vector v, with nstehfest elements, contains the coefficients
used by the Stehfest algorithm.

References:
Stehfest (1970a), Algorithm 368: Numerical inversion of Laplace transforms
[D5], Communications of the ACM, Volume 13 Issue 1, Jan. 1970, Pages 47-49

Stehfest (1970b), Remark on algorithm 368 [D5]: Numerical inversion of
Laplace transforms, Communications of the ACM, Volume 13 Issue 10, Oct. 1970,
Page 624
*/
/*---------------------------------------------------------------------------*/
double stehfest_ilt(double const (*f)(const void *, double), const void *par,
					int nstehfest, const double *const v, double t)
{
	int i;
	double
		u    = 0.0,
		sum  = 0.0,
		ln2t = log(2.0) / t;

	for(i = 0; i < nstehfest; i++){
		u = ((double)i + 1.0) * ln2t;
		sum += v[i] * f(par, u);
	}

	return (ln2t * sum);
}
