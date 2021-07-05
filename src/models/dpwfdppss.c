#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "../stehfest/stehfest.h"

#include "dpwfdppss.h"

/******** DOUBLE POROSITY INFINITE RESERVOIR, PSS INTERPOROSITY FLOW *********/
/* References:
*
* J.E. Warren and P.J. Root, The Behavior of Naturally Fractured Reservoirs,
* SPE-426-PA, 1963
*/


/**
* interporosity flow function
*/
double fpss(const double uD, const double omega, const double lambdas)
{
    double fpss;

    fpss = (uD * omega * (1.0 - omega) + lambdas)
         / (uD *         (1.0 - omega) + lambdas);

    return fpss;
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwfdppss_dCbar(const void *parameters, double u)
{
    double bc_a, dpwfb;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    bc_a  = (p->C3) / (p->C1 * p->C2 * p->qB);

    dpwfb = dpwfdpbar(p, u);

    return -bc_a * u * u * (dpwfb * dpwfb);
}
/*****************************************************************************/


/**
* deltapwf (pressure drop per unit constant flow rate) function in physical
* space. numerical inversion of the Laplace transform using Stehfest's
* algorithm
*/
double dpwfdppss(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&dpwfdpbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfdppss_dC(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwfdppss_dCbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
