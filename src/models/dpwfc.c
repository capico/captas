#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "../stehfest/stehfest.h"

#include "dpwfc.h"

/****************************** CHANNEL FLOW *********************************/
/* References:
*
* R. Nutakki and L. Mattar, Pressure Transient Analysis of Wells in Very Long
* Narrow Reservoirs,  Society of Petroleum Engineers Paper 11221, 1982
*
* A.F. Van Everdingen and W. Hurst, The Application of the Laplace
* Transformation to Flow Problems in Reservoirs, SPE-949305-G, 1949
*/


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for a homogeneous reservoir with two sealing faults, with skin factor
* and wellbore storage effects. Skin factor represented using the equivalent
* radius rws.
*/
double dpwfcbar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0  = u*b;
    aux1  = sqrt(aux0);
    aux2  = sqrt(u*z);

    numer = gsl_sf_bessel_K0(aux1);

    while(fabs(delta) > epsilon)
    {
        wn = 2.0*(p->w1*i + p->w2*j);
        ws = 2.0*(p->w1*j + p->w2*i);

        delta = gsl_sf_bessel_K0(wn * aux2)
                + gsl_sf_bessel_K0(ws * aux2);

        sum += delta;

        if(nloop%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }

        nloop++;
    }

    //printf("i = %d, j = %d\n", i, j);

    numer += sum;
    numer /= aux0;
    denom  = ( 1.0 + CD * aux0 * aux0 * numer );

    return (b * a) * (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dk function in the Laplace space
*/
double ddpwfc_dkbar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, sum1 = 0.0, delta1 = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0  = u * b;
    aux1  = sqrt(aux0);
    aux2  = sqrt(u * z);
    sum   = gsl_sf_bessel_K0(aux1);

    while( (fabs(delta) > epsilon) && (fabs(delta1) > epsilon) )
    {
        wn = 2.0*(p->w1*i + p->w2*j);
        ws = 2.0*(p->w1*j + p->w2*i);

        delta  = gsl_sf_bessel_K0(wn * aux2)
                 + gsl_sf_bessel_K0(ws * aux2);

        delta1 = wn * gsl_sf_bessel_K1(wn * aux2)
                 + ws * gsl_sf_bessel_K1(ws * aux2);

        sum  += delta;

        sum1 += delta1;

        if(nloop%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }

        nloop++;
    }

    //printf("i = %d, j = %d\n", i, j);

    numer  = (a / p->k) * ( aux1*gsl_sf_bessel_K1(aux1) + aux2*sum1 - 2.0*sum );

    denom  = ( 1.0 + CD * aux0 * sum );
    denom *= denom;
    denom *= 2.0 * u;

    return (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwfc_dCbar(const void *parameters, double u)
{
    double a, b, c, CD, rws, z, numer, denom, aux0, aux1, aux2, aux3;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);
    c   = (p->C3)        / (p->phi * p->h  * p->ct * rws * rws);

    aux0 = u * b;
    aux1 = sqrt(aux0);
    aux2 = sqrt(u * z);
    aux3 = gsl_sf_bessel_K0(aux1);

    while(fabs(delta) > epsilon)
    {
        wn = 2.0*(p->w1*i + p->w2*j);
        ws = 2.0*(p->w1*j + p->w2*i);

        delta = gsl_sf_bessel_K0(wn * aux2)
                + gsl_sf_bessel_K0(ws * aux2);

        sum += delta;

        if(nloop%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }

        nloop++;
    }

    numer  = aux3 + sum;
    numer *= numer;
    numer *= - a * b * c;

    denom  = 1.0 + aux0 * CD * (aux3 + sum);
    denom *= denom;

    return (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dS function in the Laplace space
*/
double ddpwfc_dSbar(const void *parameters, double u)
{
    double a, b, rws, w, y, z, numer, denom, aux0, aux1, aux2, aux3;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, delta = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws)     / (p->k * p->C1);
    y   = (p->phi * p->mu * p->ct * p->rw * p->rw) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)                 / (p->k * p->C1);
    w   = (p->C * p->C3) / (p->phi * p->h  * p->ct * p->rw * p->rw);

    aux0 = u * b;
    aux1 = sqrt(aux0);
    aux2 = sqrt(u * z);
    aux3 = gsl_sf_bessel_K0(aux1);

    while(fabs(delta) > epsilon)
    {
        wn = 2.0*(p->w1*i + p->w2*j);
        ws = 2.0*(p->w1*j + p->w2*i);

        delta = gsl_sf_bessel_K0(wn * aux2)
                + gsl_sf_bessel_K0(ws * aux2);

        sum += delta;

        if(nloop%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }

        nloop++;
    }

    numer  = a * aux1 * gsl_sf_bessel_K1(aux1);

    denom  = 1.0 + u * w * y * (aux3 + sum);
    denom *= denom;
    denom *= u;

    return (numer / denom);
}
/*****************************************************************************/



/**
d(dpwf)/dw1 function in the Laplace space
*/
double ddpwfc_dw1bar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, sum1 = 0.0, delta = 1.0, delta1 = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0  = u*b;
    aux1  = sqrt(aux0);
    aux2  = sqrt(u*z);

    while( (fabs(delta) > epsilon) && (fabs(delta1) > epsilon) )
    {
        wn = 2.0*(p->w1*i + p->w2*j);
        ws = 2.0*(p->w1*j + p->w2*i);

        delta  = gsl_sf_bessel_K0(wn * aux2)
                 + gsl_sf_bessel_K0(ws * aux2);

        delta1 = gsl_sf_bessel_K1(wn * aux2) * i
                 + gsl_sf_bessel_K1(ws * aux2) * j;

        sum  += delta;

        sum1 += delta1;

        if(nloop%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }

        nloop++;
    }

    //printf("i = %d, j = %d\n", i, j);

    numer  = -2.0 * a * aux2 * sum1;

    denom  = 1.0 + aux0 * CD * (gsl_sf_bessel_K0(aux1) + sum);
    denom *= denom;
    denom *= u;

    return (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dw2 function in the Laplace space
*/
double ddpwfc_dw2bar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    double epsilon = DBL_EPSILON;
    double sum = 0.0, sum1 = 0.0, delta = 1.0, delta1 = 1.0, wn, ws;
    int j = 0, i = 1, nloop = 0;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0  = u*b;
    aux1  = sqrt(aux0);
    aux2  = sqrt(u*z);

    while( (fabs(delta) > epsilon) && (fabs(delta1) > epsilon) )
    {
        wn = 2.0*(p->w1*i + p->w2*j);
        ws = 2.0*(p->w1*j + p->w2*i);

        delta  = gsl_sf_bessel_K0(wn * aux2)
                 + gsl_sf_bessel_K0(ws * aux2);

        delta1 = gsl_sf_bessel_K1(wn * aux2) * j
                 + gsl_sf_bessel_K1(ws * aux2) * i;

        sum  += delta;

        sum1 += delta1;

        if(nloop%2 == 0)
        {
            j++;
        }
        else
        {
            i++;
        }

        nloop++;
    }

    //printf("i = %d, j = %d\n", i, j);

    numer  = -2.0 * a * aux2 * sum1;

    denom  = 1.0 + aux0 * CD * (gsl_sf_bessel_K0(aux1) + sum);
    denom *= denom;
    denom *= u;

    return (numer / denom);
}
/*****************************************************************************/


/**
deltapwf function in physical space. numerical inversion of the Laplace
transform using Stehfest's algorithm
*/
double dpwfc(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&dpwfcbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dk function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfc_dk(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwfc_dkbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfc_dC(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwfc_dCbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dS function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfc_dS(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwfc_dSbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dw1 function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfc_dw1(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwfc_dw1bar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dw2 function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwfc_dw2(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwfc_dw2bar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
