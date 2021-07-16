#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "../stehfest/stehfest.h"

#include "dpwf.h"

/******************************* RADIAL FLOW *********************************/
/* References:
*
* R. Agarwal et al, An Investigation of Wellbore Storage and Skin Effect in
* Unsteady Liquid Flow: I. Analytical Treatment, SPE-2466-PA,1970
*
* A.F. Van Everdingen and W. Hurst, The Application of the Laplace
* Transformation to Flow Problems in Reservoirs, SPE-949305-G, 1949
*/


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for a infinite homogeneous reservoir with skin factor and wellbore
* storage effects. Skin factor represented using the equivalent radius
* rws = rw*exp(-S).
*/
double dpwfbar(const void *parameters, double u)
{
	double k, S, C, rws, a, b, CD, uD, uDsqrt, f;
	modelparameters *p;

	p = (modelparameters *)parameters;

	gsl_set_error_handler_off();

	k   = p->rpval[PERMEABILITY];
	S   = p->rpval[SKIN_FACTOR];
	C   = p->rpval[WELLBORE_STORAGE];

	rws = p->rw*exp(-S);
	a   = (p->qB * p->mu * p->C2)/ (p->h * k);
	b   = (p->phi * p->mu * p->ct * rws * rws) / (k * p->C1);
	CD  = (C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);
	uD  = (u * b);

	uDsqrt = sqrt(uD);
	f      = gsl_sf_bessel_K0(uDsqrt);
	f     /= uD * uDsqrt * gsl_sf_bessel_K1(uDsqrt);

	return a * b * f / (1.0 + CD * uD * uD * f);
}
/*****************************************************************************/


/**
d(dpwf)/dk function in the Laplace space
*/
double ddpwf_dkbar(const void *parameters, double u)
{
	double k, S, C, rws, a, b, y, CD, uD, uDsqrt, aux1, aux2, aux3, aux4;

	modelparameters *p;

	p = (modelparameters *)parameters;

	gsl_set_error_handler_off();

	k   = p->rpval[PERMEABILITY];
	S   = p->rpval[SKIN_FACTOR];
	C   = p->rpval[WELLBORE_STORAGE];

	rws = p->rw*exp(-S);
	a   = (p->qB * p->mu * p->C2) / (p->h * k);
	y   = (p->phi * p->mu * p->ct * rws * rws)/(p->C1);
	b   = y / k;
	CD  = (C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);
	uD  = (u * b);

	uDsqrt = sqrt(uD);
	aux2   = gsl_sf_bessel_K0(uDsqrt);
	aux3   = uDsqrt * gsl_sf_bessel_K1(uDsqrt);
	aux4   = (CD * u * b * aux2 + aux3);
	aux4   = aux4 * aux4;

	return (-a/(2.0*u* k *aux4)) * (u*b*aux2*aux2 + 2.0*aux2*aux3 - aux3*aux3);
}
/*****************************************************************************/


/**
d(dpwf)/dS function in the Laplace space
*/
double ddpwf_dSbar(const void *parameters, double u)
{
	double k, S, C, rws, uD, uDsqrt, a, b, y, z, aux1, aux2, aux3, aux4;
	modelparameters *p;

	p = (modelparameters *)parameters;

	gsl_set_error_handler_off();

	k   = p->rpval[PERMEABILITY];
	S   = p->rpval[SKIN_FACTOR];
	C   = p->rpval[WELLBORE_STORAGE];

	rws = p->rw*exp(-S);
	a   = (p->qB * p->mu * p->C2) / (p->h * k);
	y   = (p->phi * p->mu * p->ct * p->rw * p->rw) / (k * p->C1);
	b   = (p->phi * p->mu * p->ct *   rws *   rws) / (k * p->C1);
	z   = (C * p->C3) / (p->phi * p->h  * p->ct * p->rw * p->rw);
	uD  = (u * b);

	uDsqrt = sqrt(uD);
	aux2   = gsl_sf_bessel_K0(uDsqrt);
	aux3   = gsl_sf_bessel_K1(uDsqrt);
	aux4   = u * y * z * aux2 + uDsqrt * aux3;
	aux4   = aux4 * aux4;

	return a * b * (aux3*aux3 - aux2*aux2) / aux4;
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwf_dCbar(const void *parameters, double u)
{
	double bc_a, dpwfb;
	modelparameters *p;

	p = (modelparameters *)parameters;

	gsl_set_error_handler_off();

	bc_a  = (p->C3) / (p->C1 * p->C2 * p->qB);

	dpwfb = dpwfbar(p, u);

	return -bc_a * u * u * (dpwfb * dpwfb);
}
/*****************************************************************************/


/**
* deltapwf (pressure drop per unit constant flow rate) function in physical
* space. numerical inversion of the Laplace transform using Stehfest's
* algorithm
*/
double dpwf(const modelparameters *p, double t)
{
	double f;

	if(t == 0.0){
		f = 0.0;
	}
	else{
		f = stehfest_ilt(&dpwfbar, p, p->nstehfest, p->v, t);
	}

	return f;
}
/*****************************************************************************/


/**
d(dpwf)/dk function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwf_dk(const modelparameters *p, double t)
{
	double f;

	if(t == 0.0){
		f = 0.0;
	}
	else{
		f = stehfest_ilt(&ddpwf_dkbar, p, p->nstehfest, p->v, t);
	}

	return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwf_dC(const modelparameters *p, double t)
{
	double f;

	if(t == 0.0){
		f = 0.0;
	}
	else{
		f = stehfest_ilt(&ddpwf_dCbar, p, p->nstehfest, p->v, t);
	}

	return f;
}
/*****************************************************************************/


/**
d(dpwf)/dS function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwf_dS(const modelparameters *p, double t)
{
	double f;

	if(t == 0.0){
		f = 0.0;
	}
	else{
		f = stehfest_ilt(&ddpwf_dSbar, p, p->nstehfest, p->v, t);
	}

	return f;
}
/*****************************************************************************/


/**
* dr/dpi function in physical space
*/
double dr_dpi(const modelparameters *p, double t)
{
	return -1.0;
}
/*****************************************************************************/
