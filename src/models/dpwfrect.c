#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include "../stehfest/stehfest.h"

#include "dpwfrect.h"

/****************************** RECTANGULAR RESERVOIR ************************/
/* References:
*
* R. C.Earlougher Jr., Advances in Well Test Analysis, SPE Monograph 5, 1977
*
* A.F. Van Everdingen and W. Hurst, The Application of the Laplace
* Transformation to Flow Problems in Reservoirs, SPE-949305-G, 1949
*/


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for a rectangular homogeneous reservoir with skin factor and wellbore
* storage effects. Skin factor represented using the equivalent radius rws.
*/
double dpwfrectbar(const void *parameters, double u)
{
    double a, b, c, CD, rws, numer, denom, aux0;
    double k, S, C, w1x, w2x, w1y, w2y;
    double epsilon = DBL_EPSILON;
    double sum, deltax, deltay, deltaxy;
    double we, ww, ws, wn, wse_sq, wne_sq, wsw_sq, wnw_sq;
    int i, j, m, l, nloop1, nloop2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    k   = p->rpval[PERMEABILITY];
    S   = p->rpval[SKIN_FACTOR];
    C   = p->rpval[WELLBORE_STORAGE];
    w1x = p->rpval[DISTANCE_TO_FAULT_1_X];
    w2x = p->rpval[DISTANCE_TO_FAULT_2_X];
    w1y = p->rpval[DISTANCE_TO_FAULT_1_Y];
    w2y = p->rpval[DISTANCE_TO_FAULT_2_Y];

    rws = p->rw*exp(-S);
    a   = (p->qB * p->mu * p->C2) / (p->h * k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (k * p->C1);
    c   = (p->phi * p->mu * p->ct) / (k * p->C1);
    CD  = (C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0  = (u*b);
    numer = gsl_sf_bessel_K0(sqrt(aux0));

    sum   = 0.0;

    i = 1;
    j = 0;
    nloop1 = 0;
    do {
        we = 2.0*(w1x*i + w2x*j);
        ww = 2.0*(w1x*j + w2x*i);
        ws = 2.0*(w1y*i + w2y*j);
        wn = 2.0*(w1y*j + w2y*i);

        deltax = gsl_sf_bessel_K0(we * sqrt(u*c))
                 + gsl_sf_bessel_K0(ww * sqrt(u*c));

        deltay  = gsl_sf_bessel_K0(ws * sqrt(u*c))
                  + gsl_sf_bessel_K0(wn * sqrt(u*c));

        sum += (deltax + deltay);

        m = 1;
        l = 0;
        nloop2  = 0;
        do {
            ws = 2.0*(w1y*m + w2y*l);
            wn = 2.0*(w1y*l + w2y*m);

            wse_sq = we*we + ws*ws;
            wne_sq = we*we + wn*wn;
            wsw_sq = ww*ww + ws*ws;
            wnw_sq = ww*ww + wn*wn;

            deltaxy = gsl_sf_bessel_K0(sqrt(u*wse_sq*c))
                      + gsl_sf_bessel_K0(sqrt(u*wne_sq*c))
                      + gsl_sf_bessel_K0(sqrt(u*wsw_sq*c))
                      + gsl_sf_bessel_K0(sqrt(u*wnw_sq*c));

            sum += deltaxy;

            if(nloop2%2 == 0) {
                l++;
            }
            else {
                m++;
            }
            nloop2++;
        } while(fabs(deltaxy) > epsilon);


        if(nloop1%2 == 0) {
            j++;
        }
        else {
            i++;
        }
        nloop1++;
    } while(fabs(deltax) > epsilon ||
            fabs(deltay) > epsilon);

    numer += sum;
    numer /= aux0;
    denom  = ( 1.0 + CD * aux0 * aux0 * numer );

    return (a * b) * (numer / denom);
}
/*****************************************************************************/


/**
deltapwf function in physical space. numerical inversion of the Laplace
transform using Stehfest's algorithm
*/
double dpwfrect(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0) {
        f = 0.0;
    }
    else {
        f = stehfest_ilt(&dpwfrectbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
