#include <float.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include "../stehfest/stehfest.h"

#include "dpwfle.h"

/******************* INFINITE RESERVOIR, LIMITED ENTRY WELL ******************/
/* References:
*
* E. Ozkan and R. Raghavan, New Solutions for Well-Test Analysis Problems:
* Part 1 - Analytical Considerations (SPE-18615-PA) & Part 2- Computational
* Considerations and Applications (SPE-18616-PA), 1991.
*
* Kuchuk F.J. and Kirwan P.A., New Skin and Wellbore Storage Type Curves for
* Partially Penetrating Wells, SPE-11676-PA, 1987.
*/


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for an infinite conductivity well with limited entry in an infinite
* homogeneous anisotropic reservoir with skin factor and wellbore storage
* effects.
*/
double dpwflebar(const void *parameters, double u)
{
	double a, d, f, en, rws, uD, CD, hD, hwD, zD, zwD;
	double k, S, C, krat, b, zw;
	double sum, delta;
	double epsilon = DBL_EPSILON;
	int n;
	modelparameters *p;

	p = (modelparameters *)parameters;

	gsl_set_error_handler_off();

	k    = p->rpval[PERMEABILITY];
	S    = p->rpval[SKIN_FACTOR];
	C    = p->rpval[WELLBORE_STORAGE];
	krat = p->rpval[PERMEABILITY_RATIO];
	b    = p->rpval[PENETRATION_RATIO];
	zw   = p->rpval[MIDPOINT_ELEVATION];

	rws  = p->rw*exp(-S);
	hD   = (p->h  /rws) * sqrt(krat);
	hwD  = (p->h*b/rws) * sqrt(krat);
	// correlation for zD*, Horne (1995)
	zD   = 0.9069 - 0.05499*log(hwD) + 0.003745*log(hwD)*log(hwD);
	zwD  = (zw/rws) * sqrt(krat);

	a   = (p->qB * p->mu * p->C2)/ (p->h * k);
	d   = (p->phi * p->mu * p->ct * rws * rws) / (k * p->C1);
	CD  = (C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);
	uD  = u * d;

	n     = 1;
	sum   = 0.0;
	delta = 1.0;
	do{
		en     = uD + (n * n * M_PI * M_PI)/(hD * hD);
		delta  = gsl_sf_bessel_K0( sqrt(en) ) / n;
		delta *= cos(n * M_PI * (zD/hD)) * cos(n * M_PI * (zwD/hD));
		delta *= sin(n * M_PI * (0.5 * hwD/hD));

		sum   += delta;

		n++;
	} while(delta >= epsilon);

	sum *= (4.0*hD / (M_PI*hwD));

	f    = (1.0/uD) * (gsl_sf_bessel_K0(sqrt(uD)) + sum);

	return a * d * f / (1.0 + CD*uD*uD*f);
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwfle_dCbar(const void *parameters, double u)
{
	double dc_a, dpwfb;
	modelparameters *p;

	p = (modelparameters *)parameters;

	gsl_set_error_handler_off();

	dc_a  = (p->C3) / (p->C1 * p->C2 * p->qB);

	dpwfb = dpwflebar(p, u);

	return -dc_a * u * u * (dpwfb * dpwfb);
}
/*****************************************************************************/


/**
* deltapwf (pressure drop per unit constant flow rate) function in physical
* space. numerical inversion of the Laplace transform using Stehfest's
* algorithm
*/
double dpwfle(const modelparameters *p, double t)
{
	double f;

	if(t == 0.0) {
		f = 0.0;
	}
	else {
		f = stehfest_ilt(&dpwflebar, p, p->nstehfest, p->v, t);
	}

	return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfle_dC(const modelparameters *p, double t)
{
	double f;

	if(t == 0.0) {
		f = 0.0;
	}
	else {
		f = stehfest_ilt(&ddpwfle_dCbar, p, p->nstehfest, p->v, t);
	}

	return f;
}
/*****************************************************************************/
