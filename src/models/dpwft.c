#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>

#include "dpwft.h"

/***************** RADIAL FLOW, TRANSFORMED PARAMETERS ***********************/
/* References:
*
* A. Dastan and R. N. Horne, Robust Well Test Interpretation by Nonlinear
* Regression with Parameter and Data Transformations, SPE Western Regional
* Meeting, SPE paper 132467, 2010
*/

/*
minimization over the transformed parameters:
kt  = ln(k),
Ct  = ln(C),
pit = ln(pi),
St  = ln(SS),
SS  = (S + 8)/(25 - S)
*/


/**
d(dpwf)/dkt function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwft_dk(const modelparameters *p, double t)
{
    return p->k * ddpwf_dk(p, t);
}
/*****************************************************************************/


/**
d(dpwf)/dCt function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwft_dC(const modelparameters *p, double t)
{
    return p->C * ddpwf_dC(p, t);
}
/*****************************************************************************/


/**
d(dpwf)/dSt function function in physical space. numerical
inversion of the Laplace transform using Stehfest's
algorithm
*/
double ddpwft_dS(const modelparameters *p, double t)
{
    return -(1.0/33.0) * (p->S - 25.0) * (8.0 + p->S) * ddpwf_dS(p, t);
}
/*****************************************************************************/


/**
* dr/dpit function in physical space
*/
double drt_dpi(const modelparameters *p, double t)
{
    return -( p->pi );
}
/*****************************************************************************/
