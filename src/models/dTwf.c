#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>

#include "dTwf.h"

/******************************* RADIAL FLOW *********************************/
/* References:
*
* Palabiyik, Y., Onur, M., Tureyen, O. I., & Cinar, M. (2016). 
* Transient temperature behavior and analysis of single-phase liquid-water 
* geothermal reservoirs during drawdown and buildup tests: Part I. 
* Theory, new analytical and approximate solutions. 
* Journal of Petroleum Science and Engineering, 146, 637–656. 
* https://doi.org/10.1016/j.petrol.2016.08.003
*/


/**
* delta Twf (temeprature drop per unit constant flow rate) function 
* for a infinite homogeneous reservoir with skin factor. 
8 Skin factor represented using the equivalent radius
* rws = rw*exp(-S).
*/
double dTwfcl(const modelparameters *p, double t)
{
    double a, b, f, d;
    double k, S, cpt;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    cpt = p->rpval[EFFECTIVE_HEAT_CAPACITY];
    a   = (p->C2 * p->qB * p->mu) / (k * p->h);
    b   = (pow(p->rw, 2.0) * p->phi * p->ct * rp->mu * p->ct) / (4.0 * k * p->C1);
    f   = (p->phi * p->rhosc / p->B * p->cp) / cpt * (p->ejt + 1 / (p->rhosc / p->B * p->cp)) - p->ejt;
    d   = (p->phi * p->rhosc / p->B * p->cp * p->ct) / p->cp / 2 * a;

    return - a * 0.5 * (-p->ejt * (gsl_sf_expint_E1(b / t) + 2 * S) - f * gsl_sf_expint_E1(b / t - d));
}
/*****************************************************************************/