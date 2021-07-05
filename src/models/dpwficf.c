#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "../stehfest/stehfest.h"

#include "dpwficf.h"

/************ INFINITE RESERVOIR INFINITE CONDUCTIVITY FRACTURE **************/
/* References:
*
* E. Ozkan and R. Raghavan, New Solutions for Well-Test Analysis Problems:
* Part 1 - Analytical Considerations (SPE-18615-PA) & Part 2- Computational
* Considerations and Applications (SPE-18616-PA), 1991.
*/


double integrand(double x, void * p)
{
    return gsl_sf_bessel_K0(x) + log(x);
}


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for an infinite conductivity fractured well in an infinite homogeneous
* reservoir with skin factor and wellbore storage effects.
*/
double dpwficfbar(const void *parameters, double u)
{
    double a, b, f, xfs, uD, CD, xD, integral_1, integral_2, lim_1, lim_2;
    double
        epsabs =  1.0e-10, // absolute error limit
        epsrel = 1.0e-6;   // relative error limit
    double abserr; // estimate of the absolute error
    size_t neval;  // number of function evaluations
    gsl_function F;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    xD  = 0.732; // for infinite conductivity behavior

    xfs = p->xf*exp(-p->S);
    a   = (p->qB * p->mu * p->C2)/ (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * xfs * xfs) / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * xfs * xfs);
    uD  = u * b;

    lim_1 = sqrt(uD) * (1.0 - xD);
    lim_2 = sqrt(uD) * (1.0 + xD);

    integral_1  = 0.0;
    integral_2  = 0.0;

    // to avoid the singularity at zero, the original integrand is written as
    // K0(x) = (K0(x) + log(x)) - log(x). the second part is analytically
    // integrated.
    F.function = &integrand;

    gsl_integration_qng(&F, 0.0, lim_1, epsabs, epsrel, &integral_1, &abserr, &neval);
    gsl_integration_qng(&F, 0.0, lim_2, epsabs, epsrel, &integral_2, &abserr, &neval);

    //printf("%E\n", abserr);
    //printf("%zu\n", neval);

    integral_1 += sqrt(uD) * (1.0 - xD) * ( 1.0 - log( sqrt(uD)*(1.0 - xD) ) );
    integral_2 += sqrt(uD) * (1.0 + xD) * ( 1.0 - log( sqrt(uD)*(1.0 + xD) ) );

    f = (integral_1 + integral_2) / (2.0 * uD * sqrt(uD)) ;

    return a * b * f / (1.0 + CD*uD*uD*f);
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwficf_dCbar(const void *parameters, double u)
{
    double bc_a, dpwfb;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    bc_a  = (p->C3) / (p->C1 * p->C2 * p->qB);

    dpwfb = dpwficfbar(p, u);

    return -bc_a * u * u * (dpwfb * dpwfb);
}
/*****************************************************************************/


/**
* deltapwf (pressure drop per unit constant flow rate) function in physical
* space. numerical inversion of the Laplace transform using Stehfest's
* algorithm
*/
double dpwficf(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&dpwficfbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwficf_dC(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwficf_dCbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
