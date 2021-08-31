#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>

#include "dTwf.h"

/******************* TEMPERATURE, DRAWDOWN, RADIAL FLOW **********************/
/* References:
*
* Palabiyik, Y., Onur, M., Tureyen, O. I., & Cinar, M. (2016).
* Transient temperature behavior and analysis of single-phase liquid-water
* geothermal reservoirs during drawdown and buildup tests: Part I.
* Theory, new analytical and approximate solutions.
* Journal of Petroleum Science and Engineering, 146, 637�656.
* https://doi.org/10.1016/j.petrol.2016.08.003
*/


/**
* delta Twf (temperature drop per unit constant flow rate) function
* for a infinite homogeneous reservoir with skin factor.
*/
double dTwf(const modelparameters *p, double t)
{
    double a, b, cpR, f, d;
    double k, S, cpt, ejt;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    cpt = p->rpval[EFFECTIVE_HEAT_CAPACITY];
    ejt = p->rpval[COEFFICIENT_JOULE_THOMSON];
    phi = p->rpval[POROSITY];

    a   = (p->C2 * p->qB * p->mu) / (k * p->h);
    b   = (p->rw * p->rw * phi * p->ct * p->mu ) / (4.0 * k * p->C1);
    cpR = (p->rhosc / p->B) * p->cp / cpt;
    f   = phi * cpR * (ejt + 1.0 / (p->rhosc * p->cp / p->B)) - ejt;
    d   = (phi * cpR * p->ct)  * a * 0.5;

    return -0.5*a*( -ejt*(gsl_sf_expint_E1(b/t) + 2.0*S) - f*gsl_sf_expint_E1(b/t + d) );

}
/*****************************************************************************/


/**
*/
double dr_dTi(const modelparameters *p, double t)
{
	return -1.0;
}
/*****************************************************************************/
