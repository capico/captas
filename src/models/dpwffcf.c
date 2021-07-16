#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "../stehfest/stehfest.h"

#include "dpwffcf.h"

/************** INFINITE RESERVOIR FINITE CONDUCTIVITY FRACTURE **************/
/* References:
*
* Lee, S. T. and Brockenbrough, I.R., A new Approximate Analytic Solution for
* Finite Conductivity Vertical Fractures, SPE-12013-PA (1986)
*
* T.A. Blasingame and B.D. Poe, Semianalytic Solutions for a Well With a Single
* Finite-Conductivity Vertical Fracture, SPE-26424-MS, 1993
*/


double integrand_inf(double x, void * p)
{
    return gsl_sf_bessel_K0(x) + log(x);
}


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for an finite conductivity fractured well in an infinite homogeneous
* reservoir with skin factor and wellbore storage effects. "desuperposition" of
* the trilinear model plus the infinite conductivity model.
*/
double dpwffcfbar(const void *parameters, double u)
{
    double a, b, f, uD, CD, xD, integral_1, integral_2, lim_1, lim_2;
    double k, S, C, xf, fc;
    double CfD, CfD_inf, psi, psi_inf;
    double
    epsabs = 1.0e-12, // absolute error limit
    epsrel = 1.0e-8;  // relative error limit
    double abserr; // estimate of the absolute error
    size_t neval;  // number of function evaluations
    gsl_function F;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
	S   = p->rpval[SKIN_FACTOR];
	C   = p->rpval[WELLBORE_STORAGE];
	xf  = p->rpval[FRACTURE_HALF_LENGTH];
	fc  = p->rpval[FRACTURE_CONDUCTIVITY];

    a   = (p->qB * p->mu * p->C2)/ (p->h * k);
    b   = (p->phi * p->mu * p->ct * xf * xf) / (k * p->C1);
    CD  = (C * p->C3) / (p->phi * p->h  * p->ct * xf * xf);
    uD  = u * b;

    CfD     = fc / (k * xf);
    CfD_inf = M_PI * 1.0e6;

    psi     = sqrt( 2.0 * (uD + sqrt(uD))/CfD     );
    psi_inf = sqrt( 2.0 * (uD + sqrt(uD))/CfD_inf );

    f  = (M_PI / uD );
    f *= (1.0 / (CfD*psi*tanh(psi)) - 1.0 / (CfD_inf*psi_inf*tanh(psi_inf)));

    xD  = 0.732; // for infinite conductivity behavior

    lim_1 = sqrt(uD) * (1.0 - xD);
    lim_2 = sqrt(uD) * (1.0 + xD);

    integral_1  = 0.0;
    integral_2  = 0.0;

    // to avoid the singularity at zero, the original integrand is written as
    // K0(x) = (K0(x) + log(x)) - log(x). the second part is analytically
    // integrated.
    F.function = &integrand_inf;

    gsl_integration_qng(&F, 0.0, lim_1, epsabs, epsrel, &integral_1, &abserr, &neval);
    gsl_integration_qng(&F, 0.0, lim_2, epsabs, epsrel, &integral_2, &abserr, &neval);

    //printf("%E\n", abserr);
    //printf("%zu\n", neval);

    integral_1 += sqrt(uD) * (1.0 - xD) * ( 1.0 - log( sqrt(uD)*(1.0 - xD) ) );
    integral_2 += sqrt(uD) * (1.0 + xD) * ( 1.0 - log( sqrt(uD)*(1.0 + xD) ) );

    f += (integral_1 + integral_2) / (2.0 * uD * sqrt(uD)) + S/uD;

    return a * b * f / (1.0 + CD*uD*uD*f);
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwffcf_dCbar(const void *parameters, double u)
{
    double bc_a, dpwfb;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    bc_a  = (p->C3) / (p->C1 * p->C2 * p->qB);

    dpwfb = dpwffcfbar(p, u);

    return -bc_a * u * u * (dpwfb * dpwfb);
}
/*****************************************************************************/


/**
d(dpwf)/dS function in the Laplace space
*/
double ddpwffcf_dSbar(const void *parameters, double u)
{
    double a, b, f, uD, CD, xD, integral_1, integral_2, lim_1, lim_2;
    double k, S, C, xf, fc;
    double CfD, CfD_inf, psi, psi_inf;
    double
    epsabs = 1.0e-12, // absolute error limit
    epsrel = 1.0e-8;  // relative error limit
    double abserr; // estimate of the absolute error
    size_t neval;  // number of function evaluations
    gsl_function F;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
	S   = p->rpval[SKIN_FACTOR];
	C   = p->rpval[WELLBORE_STORAGE];
	xf  = p->rpval[FRACTURE_HALF_LENGTH];
	fc  = p->rpval[FRACTURE_CONDUCTIVITY];

    a   = (p->qB * p->mu * p->C2)/ (p->h * k);
    b   = (p->phi * p->mu * p->ct * xf * xf) / (k * p->C1);
    CD  = (C * p->C3) / (p->phi * p->h  * p->ct * xf * xf);
    uD  = u * b;

    CfD     = fc / (k * xf);
    CfD_inf = M_PI * 1.0e6;

    psi     = sqrt( 2.0 * (uD + sqrt(uD))/CfD     );
    psi_inf = sqrt( 2.0 * (uD + sqrt(uD))/CfD_inf );

    f  = (M_PI / uD );
    f *= (1.0 / (CfD*psi*tanh(psi)) - 1.0 / (CfD_inf*psi_inf*tanh(psi_inf)));

    xD  = 0.732; // for infinite conductivity behavior

    lim_1 = sqrt(uD) * (1.0 - xD);
    lim_2 = sqrt(uD) * (1.0 + xD);

    integral_1  = 0.0;
    integral_2  = 0.0;

    // to avoid the singularity at zero, the original integrand is written as
    // K0(x) = (K0(x) + log(x)) - log(x). the second part is analytically
    // integrated.
    F.function = &integrand_inf;

    gsl_integration_qng(&F, 0.0, lim_1, epsabs, epsrel, &integral_1, &abserr, &neval);
    gsl_integration_qng(&F, 0.0, lim_2, epsabs, epsrel, &integral_2, &abserr, &neval);

    integral_1 += sqrt(uD) * (1.0 - xD) * ( 1.0 - log( sqrt(uD)*(1.0 - xD) ) );
    integral_2 += sqrt(uD) * (1.0 + xD) * ( 1.0 - log( sqrt(uD)*(1.0 + xD) ) );

    f += (integral_1 + integral_2) / (2.0 * uD * sqrt(uD)) + S/uD;

    return a * b / ( uD * (1.0 + CD*uD*uD*f) * (1.0 + CD*uD*uD*f) );
}
/*****************************************************************************/


/**
* deltapwf (pressure drop per unit constant flow rate) function in physical
* space. numerical inversion of the Laplace transform using Stehfest's
* algorithm
*/
double dpwffcf(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0) {
        f = 0.0;
    }
    else {
        f = stehfest_ilt(&dpwffcfbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwffcf_dC(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0) {
        f = 0.0;
    }
    else {
        f = stehfest_ilt(&ddpwffcf_dCbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dS function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwffcf_dS(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0) {
        f = 0.0;
    }
    else {
        f = stehfest_ilt(&ddpwffcf_dSbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
