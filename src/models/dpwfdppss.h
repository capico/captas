#ifndef DPWFDPPSS_H
#define DPWFDPPSS_H

#include "modelparam.h"
#include "dpwfdp.h"

double fpss(const double uD, const double omega, const double lambdas);

double ddpwfdppss_dCbar(const void *parameters, double u);

double dpwfdppss(const modelparameters *p, double t);

double ddpwfdppss_dC(const modelparameters *p, double t);

#endif
