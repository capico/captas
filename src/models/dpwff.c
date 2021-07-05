#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include "../stehfest/stehfest.h"

#include "dpwff.h"

/****************************** SEALING FAULT ********************************/
/* References:
*
* R. C.Earlougher Jr., Advances in Well Test Analysis, SPE Monograph 5, 1977
*
* A.F. Van Everdingen and W. Hurst, The Application of the Laplace
* Transformation to Flow Problems in Reservoirs, SPE-949305-G, 1949
*/


/**
delta pwf function in the Laplace space, for a homogeneous reservoir, with
skin factor and wellbore storage effects, and the well close to a linear
sealing fault.
*/
double dpwffbar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0 = (u * b);
    aux1 = gsl_sf_bessel_K0(sqrt(aux0));
    aux2 = gsl_sf_bessel_K0(2.0 * p->L * sqrt(u * z));

    numer = (aux1 + aux2) / aux0;

    denom = ( 1.0 + aux0 * aux0 * CD * numer );

    return (a * b) * (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dk function in the Laplace space
*/
double ddpwff_dkbar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0 = (u * b);
    aux1 = gsl_sf_bessel_K0(sqrt(aux0));
    aux2 = gsl_sf_bessel_K0(2.0 * p->L * sqrt(u * z));

    numer = (a / p->k) * ( sqrt(aux0) * gsl_sf_bessel_K1(sqrt(aux0))
                           + 2.0 *p->L * sqrt(u * z) * gsl_sf_bessel_K1(2.0 * p->L * sqrt(u * z))
                           - 2.0 * (aux1 + aux2) );

    denom  = ( 1.0 + aux0 * CD * (aux1 + aux2) );
    denom *= denom;
    denom  = 2.0 * u * denom;


    return (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dC function in the Laplace space
*/
double ddpwff_dCbar(const void *parameters, double u)
{
    double a, b, c, rws, z, numer, denom, aux0, aux1, aux2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    c   = p->C3 / (p->phi * p->h  * p->ct * rws * rws);

    aux0 = (u * b);
    aux1 = gsl_sf_bessel_K0(sqrt(aux0));
    aux2 = gsl_sf_bessel_K0(2.0 * p->L * sqrt(u * z));

    numer  = (aux1 + aux2);
    numer *= numer;
    numer *= -a * b * c;

    denom  = ( 1.0 + b * u * c * p->C * (aux1 + aux2) );
    denom *= denom;

    return (numer / denom);
}
/*****************************************************************************/


/**
d(dpwf)/dS function in the Laplace space
*/
double ddpwff_dSbar(const void *parameters, double u)
{
    double a, b, rws, w, y, z, numer, denom, aux0, aux1, aux2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws)     / (p->k * p->C1);
    y   = (p->phi * p->mu * p->ct * p->rw * p->rw) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)                 / (p->k * p->C1);
    w   = (p->C * p->C3) / (p->phi * p->h  * p->ct * p->rw * p->rw);

    aux0 = (u * b);
    aux1 = gsl_sf_bessel_K0(sqrt(aux0));
    aux2 = gsl_sf_bessel_K0(2.0 * p->L * sqrt(u * z));

    numer  = a * sqrt(aux0) * gsl_sf_bessel_K1(sqrt(aux0));

    denom  = ( 1.0 + u * y * w * (aux1 + aux2) );
    denom *= denom;
    denom *= u;

    return (numer / denom);

}
/*****************************************************************************/


/**
d(dpwf)/dL function in the Laplace space
*/
double ddpwff_dLbar(const void *parameters, double u)
{
    double a, b, CD, rws, z, numer, denom, aux0, aux1, aux2;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = p->rw*exp(-p->S);
    a   = (p->qB * p->mu * p->C2) / (p->h * p->k);
    b   = (p->phi * p->mu * p->ct * rws * rws) / (p->k * p->C1);
    z   = (p->phi * p->mu * p->ct)             / (p->k * p->C1);
    CD  = (p->C * p->C3) / (p->phi * p->h  * p->ct * rws * rws);

    aux0 = (u * b);
    aux1 = gsl_sf_bessel_K0(sqrt(aux0));
    aux2 = gsl_sf_bessel_K0(2.0 * p->L * sqrt(u * z));

    denom  = 1.0 + aux0 * CD * (aux1 + aux2);
    denom *= denom;
    denom *= u;

    numer  = -2.0 * a * sqrt(u * z) * gsl_sf_bessel_K1(2.0 * p->L * sqrt(u * z));

    return (numer / denom);
}
/*****************************************************************************/


/**
deltapwf function in physical space. numerical inversion of the Laplace
transform using Stehfest's algorithm
*/
double dpwff(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&dpwffbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dk function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwff_dk(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwff_dkbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwff_dC(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwff_dCbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dS function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwff_dS(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwff_dSbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dpwf)/dL function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwff_dL(const modelparameters *p, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddpwff_dLbar, p, p->nstehfest, p->v, t);
    }

    return f;
}
/*****************************************************************************/
