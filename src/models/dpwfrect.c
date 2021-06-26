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
    double epsilon = DBL_EPSILON;
    double sum, deltax, deltay, deltaxy;
    double we, ww, ws, wn, wse_sq, wne_sq, wsw_sq, wnw_sq;
    int i, j, k, l, nloop1, nloop2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    c   = (p->phi * p->mu * p->ct) / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0  = (u*b);
    numer = gsl_sf_bessel_K0(sqrt(aux0));

    sum   = 0.0;

    i = 1;
    j = 0;
    nloop1 = 0;
    do
    {
        we = 2.0*(p->w1x*i + p->w2x*j);
        ww = 2.0*(p->w1x*j + p->w2x*i);
        ws = 2.0*(p->w1y*i + p->w2y*j);
        wn = 2.0*(p->w1y*j + p->w2y*i);

        deltax = gsl_sf_bessel_K0(we * sqrt(u*c))
                 + gsl_sf_bessel_K0(ww * sqrt(u*c));

        deltay  = gsl_sf_bessel_K0(ws * sqrt(u*c))
                  + gsl_sf_bessel_K0(wn * sqrt(u*c));

        sum += (deltax + deltay);

        k = 1;
        l = 0;
        nloop2  = 0;
        do
        {
            ws = 2.0*(p->w1y*k + p->w2y*l);
            wn = 2.0*(p->w1y*l + p->w2y*k);

            wse_sq = we*we + ws*ws;
            wne_sq = we*we + wn*wn;
            wsw_sq = ww*ww + ws*ws;
            wnw_sq = ww*ww + wn*wn;

            deltaxy = gsl_sf_bessel_K0(sqrt(u*wse_sq*c))
                      + gsl_sf_bessel_K0(sqrt(u*wne_sq*c))
                      + gsl_sf_bessel_K0(sqrt(u*wsw_sq*c))
                      + gsl_sf_bessel_K0(sqrt(u*wnw_sq*c));

            sum += deltaxy;

            if(nloop2%2 == 0)
            {
                l++;
            }
            else
            {
                k++;
            }
            nloop2++;
        }
        while(fabs(deltaxy) > epsilon);


        if(nloop1%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }
        nloop1++;
    }
    while(fabs(deltax) > epsilon ||
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

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&dpwfrectbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
