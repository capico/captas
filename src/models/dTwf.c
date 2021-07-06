#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "../stehfest/stehfest.h"

#include "dTwf.h"

/******************************* RADIAL FLOW *********************************/
/* References:
*
* Dr. Onur's paper
*/


/**
* delta Twf (pressure drop per unit constant flow rate) function in the Laplace
* space, for a infinite homogeneous reservoir with skin factor and wellbore
* storage effects. Skin factor represented using the equivalent radius
* rws = rw*exp(-S).
*/
double dTwfbar(const void *parameters, double u)
{
    double rws, a, b, CD, aux1, aux2, aux3;
    modelparameters *T;

    T = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = T->rw*exp(-T->S);
    a   = (T->qB * T->mu * T->C2)/ (T->h * T->k);
    b   = (T->phi * T->mu * T->ct * rws * rws) / (T->k * T->C1);
    CD  = (T->C * T->C3) / (T->phi * T->h  * T->ct * rws * rws);

    aux1  = sqrt(u * b);
    aux2  = gsl_sf_bessel_K0(aux1);
    aux3  = aux1 * gsl_sf_bessel_K1(aux1);

    return a*b*aux2 / (CD*u*u*b*b*aux2 +  u*b*aux3);
}
/*****************************************************************************/


/**
d(dTwf)/dk function in the Laplace space
*/
double ddTwf_dkbar(const void *parameters, double u)
{
    double rws, a, b, y, CD, aux1, aux2, aux3, aux4;

    modelparameters *T;

    T = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = T->rw*exp(-T->S);
    a   = (T->qB * T->mu * T->C2) / (T->h * T->k);
    y   = (T->phi * T->mu * T->ct * rws * rws)/(T->C1);
    b   = y / T->k;
    CD  = (T->C * T->C3) / (T->phi * T->h  * T->ct * rws * rws);

    aux1  = sqrt(u * b);
    aux2  = gsl_sf_bessel_K0(aux1);
    aux3  = aux1 * gsl_sf_bessel_K1(aux1);
    aux4  = (CD * u * b * aux2 + aux3);
    aux4  = aux4 * aux4;

    return (-a/(2.0*u*T->k*aux4)) * (u*b*aux2*aux2 + 2.0*aux2*aux3 - aux3*aux3);
}
/*****************************************************************************/


/**
d(dTwf)/dC function in the Laplace space
*/
//double ddTwf_dCbar(const void *parameters, double u)
//{
//    double rws, a, b, z, CD, aux1, aux2, aux3, aux4;
//    modelparameters *T;
//
//	T = (modelparameters *)parameters;
//
//	gsl_set_error_handler_off();
//
//	rws = T->rw*exp(-v->S);
//	a   = (T->qB * T->mu * T->C2) / (T->h * T->k);
//	b   = (T->phi * T->mu * T->ct * rws * rws)/(T->k * T->C1);
//	z   = (T->C3)        / (T->phi * T->h  * T->ct * rws * rws);
//	CD  = T->C * z;
//
//	aux1  = sqrt(u * b);
//	aux2  = gsl_sf_bessel_K0(aux1);
//	aux3  = aux1 * gsl_sf_bessel_K1(aux1);
//	aux4  = CD * u * b * aux2 + aux3;
//	aux4  = aux4 * aux4;
//
//	return - a * b * z * aux2 * aux2 / aux4;
//}
/*****************************************************************************/


/**
d(dTwf)/dS function in the Laplace space
*/
double ddTwf_dSbar(const void *parameters, double u)
{
    double rws, a, b, y, z, aux1, aux2, aux3, aux4;
    modelparameters *T;

    T = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    rws = T->rw*exp(-T->S);
    a   = (T->qB * T->mu * T->C2) / (T->h * T->k);
    y   = (T->phi * T->mu * T->ct * T->rw * T->rw) / (T->k * T->C1);
    b   = (T->phi * T->mu * T->ct * rws * rws)     / (T->k * T->C1);
    z   = (T->C * T->C3) / (T->phi * T->h  * T->ct * T->rw * T->rw);

    aux1  = sqrt(u * b);
    aux2  = gsl_sf_bessel_K0(aux1);
    aux3  = gsl_sf_bessel_K1(aux1);
    aux4  = u * y * z * aux2 + aux1 * aux3;
    aux4  = aux4 * aux4;

    return a * b * (aux3*aux3 - aux2*aux2) / aux4;
}
/*****************************************************************************/


/**
d(dTwf)/dC function in the Laplace space
*/
double ddTwf_dCbar(const void *parameters, double u)
{
    double bc_a, dTwfb;
    modelparameters *T;

    T = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    bc_a  = (T->C3) / (T->C1 * T->C2 * T->qB);

    dTwfb = dTwfbar(T, u);

    return -bc_a * u * u * (dTwfb * dTwfb);
}
/*****************************************************************************/


/**
* deltaTwf (pressure drop per unit constant flow rate) function in physical
* space. numerical inversion of the Laplace transform using Stehfest's
* algorithm
*/
double dTwf(const modelparameters *T, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&dTwfbar, T, T->nstehfest, T->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dTwf)/dk function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddTwf_dk(const modelparameters *T, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddTwf_dkbar, T, T->nstehfest, T->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dTwf)/dC function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddTwf_dC(const modelparameters *T, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddTwf_dCbar, T, T->nstehfest, T->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
d(dTwf)/dS function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddTwf_dS(const modelparameters *T, double t)
{
    double f;

    if(t == 0.0)
    {
        f = 0.0;
    }
    else
    {
        f = stehfest_ilt(&ddTwf_dSbar, T, T->nstehfest, T->v, t);
    }

    return f;
}
/*****************************************************************************/


/**
* dr/dpi function in physical space
*/
double dr_dpi(const modelparameters *T, double t)
{
    return -1.0;
}
/*****************************************************************************/
