#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>

#include "dpwfcl.h"

/*********************** CHANNEL FLOW - NO WELLBORE STORAGE ******************/
/* References:
*
* R. Nutakki and L. Mattar, Pressure Transient Analysis of Wells in Very Long
* Narrow Reservoirs,  Society of Petroleum Engineers Paper 11221, 1982
*/


/**
* delta pwf (pressure drop per unit constant flow rate) function, for a homo-
* geneous reservoir with two sealing faults, with skin factor effects.
*/
double dpwfcl(const modelparameters *p, double t)
{
    double a, b, rws, z;
    double k, S, w1, w2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    w1  = p->rpval[DISTANCE_TO_FAULT_1];
    w2  = p->rpval[DISTANCE_TO_FAULT_2];
    phi = p->rpval[POROSITY];

    rws = p->rw * gsl_sf_exp(-S);
    a   = (p->qB * p->mu * p->C2) / (p->h * k);
    b   = (phi * p->mu * p->ct * rws * rws) / (4.0*k * p->C1 * t);
    z   = (phi * p->mu * p->ct)             / (4.0*k * p->C1 * t);

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


/**
d(dpwf)/dk
*/
double ddpwfcl_dk(const modelparameters *p, double t)
{
    double a, b, rws, z;
    double k, S, w1, w2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, sum1 = 0.0, delta1 = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    w1  = p->rpval[DISTANCE_TO_FAULT_1];
    w2  = p->rpval[DISTANCE_TO_FAULT_2];
    phi = p->rpval[POROSITY];

    rws = p->rw * gsl_sf_exp(-S);
    a   = (p->qB * p->mu * p->C2) / (p->h * k);
    b   = (phi * p->mu * p->ct * rws * rws) / (4.0*k * p->C1 * t);
    z   = (phi * p->mu * p->ct)             / (4.0*k * p->C1 * t);

    while( (fabs(delta) > epsilon) && (fabs(delta1) > epsilon) ) {
        wn = 2.0*(w1*i + w2*j);
        ws = 2.0*(w1*j + w2*i);

        delta  = gsl_sf_expint_E1(wn*wn*z)
                 + gsl_sf_expint_E1(ws*ws*z);

        delta1 = gsl_sf_exp(-wn*wn*z)
                 + gsl_sf_exp(-ws*ws*z);

        sum  += delta;

        sum1 += delta1;

        if(nloop%2 == 0) {
            j++;
        }
        else {
            i++;
        }

        nloop++;
    }

    return -(a/k) * 0.5 * (gsl_sf_expint_E1(b) + sum - sum1);
}
/*****************************************************************************/


/**
d(dpwf)/dS
*/
double ddpwfcl_dS(const modelparameters *p, double t)
{
    double a, b, rws;
    double k, S;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    phi = p->rpval[POROSITY];

    rws = p->rw * gsl_sf_exp(-S);
    a   = (p->qB * p->mu * p->C2) / (p->h * k);
    b   = (phi * p->mu * p->ct * rws * rws) / (4.0*k * p->C1 * t);

    return a * gsl_sf_exp(-b);
}
/*****************************************************************************/


/**
d(dpwf)/dw1
*/
double ddpwfcl_dw1(const modelparameters *p, double t)
{
    double a, z;
    double k, w1, w2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    w1  = p->rpval[DISTANCE_TO_FAULT_1];
    w2  = p->rpval[DISTANCE_TO_FAULT_2];
    phi = p->rpval[POROSITY];

    a = (p->qB * p->mu * p->C2) / (p->h * k);
    z = (phi * p->mu * p->ct) / (4.0*k * p->C1 * t);

    while(fabs(delta) > epsilon) {
        wn = 2.0*(w1*i + w2*j);
        ws = 2.0*(w1*j + w2*i);

        delta = (gsl_sf_exp(-wn*wn*z) / wn) * i
                + (gsl_sf_exp(-ws*ws*z) / ws) * j;

        sum += delta;

        if(nloop%2 == 0) {
            j++;
        }
        else {
            i++;
        }

        nloop++;
    }

    return -a * 2.0 * sum;
}
/*****************************************************************************/


/**
d(dpwf)/dw2
*/
double ddpwfcl_dw2(const modelparameters *p, double t)
{
    double a, z;
    double k, w1, w2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    w1  = p->rpval[DISTANCE_TO_FAULT_1];
    w2  = p->rpval[DISTANCE_TO_FAULT_2];
    phi = p->rpval[POROSITY];

    a = (p->qB * p->mu * p->C2) / (p->h * k);
    z = (phi * p->mu * p->ct) / (4.0*k * p->C1 * t);

    while(fabs(delta) > epsilon) {
        wn = 2.0*(w1*i + w2*j);
        ws = 2.0*(w1*j + w2*i);

        delta = (gsl_sf_exp(-wn*wn*z) / wn) * j
                + (gsl_sf_exp(-ws*ws*z) / ws) * i;

        sum += delta;

        if(nloop%2 == 0) {
            j++;
        }
        else {
            i++;
        }

        nloop++;
    }

    return -a * 2.0 * sum;
}
/*****************************************************************************/
