#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "../stehfest/stehfest.h"

#include "dpwfdp.h"

/*********************** DOUBLE POROSITY INFINITE RESERVOIR ******************/
/* References:
*
* J.E. Warren and P.J. Root, The Behavior of Naturally Fractured Reservoirs,
* SPE-426-PA, 1963
*
* Bourdet, D. and Gringarten, A.C., Determination of Fissured Volume and Block
* Size in Fractured Reservoirs by Type-Curve Analysis, SPE 9293,1980
*/


/**
* delta pwf (pressure drop per unit constant flow rate) function in the Laplace
* space, for a double porosity infinite reservoir with skin factor and wellbore
* storage effects.
*/
double dpwfdpbar(const void *parameters, double u)
{
    double a, b, rws, lambdas, CD, f, uD, aux1, aux2;
    double k, S, C, omega, lambda;
    modelparameters *p;

    p = (modelparameters *)parameters;

    gsl_set_error_handler_off();

    k      = p->rpval[PERMEABILITY];
	S      = p->rpval[SKIN_FACTOR];
	C      = p->rpval[WELLBORE_STORAGE];
	omega  = p->rpval[OMEGA];
	lambda = p->rpval[LAMBDA];
    phi = p->rpval[POROSITY];

    rws     = p->rw*exp(-S);
    lambdas = lambda*exp(-2.0*S);

    a   = (p->qB * p->mu * p->C2)/ (p->h * k);
    b   = (phi * p->mu * p->ct * rws * rws) / (k * p->C1);
    CD  = (C * p->C3) / (phi * p->h  * p->ct * rws * rws);

    uD  = u * b;
    f   = (p->f)[p->model](uD, omega, lambdas);

    aux1  = sqrt(uD * f);
    aux2  = gsl_sf_bessel_K0(aux1);
    aux2 /= uD * aux1 * gsl_sf_bessel_K1(aux1);

    return a * b * aux2 / (1.0 +  CD*uD*uD * aux2);
}
/*****************************************************************************/

