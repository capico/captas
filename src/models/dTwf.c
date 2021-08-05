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
    double a, b, rws, z;
    double k, S, cpt;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    cpt = p->rpval[EFFECTIVE_HEAT_CAPACITY];
    rws = p->rw * gsl_sf_exp(-S);
    a   = (p->qB * p->mu * p->C2) / (p->h * k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (4.0*k * p->C1 * t);
    z   = (p->phi * p->mu * p->ct)             / (4.0*k * p->C1 * t);

    while(fabs(delta) > epsilon) {
        wn = 2.0*(w1*i + w2*j);
        ws = 2.0*(w1*j + w2*i);

        delta = gsl_sf_expint_E1(wn*wn*z)
                + gsl_sf_expint_E1(ws*ws*z);

        sum += delta;

        if(nloop%2 == 0) {
            j++;
        }
        else {
            i++;
        }

        nloop++;
    }

    return a * 0.5 * (gsl_sf_expint_E1(b) + sum);
}
/*****************************************************************************/