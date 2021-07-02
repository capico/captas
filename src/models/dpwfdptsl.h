#ifndef DPWFDPTSL_H
#define DPWFDPTSL_H

#include "modelparam.h"
#include "dpwfdp.h"

double ftsl(const double uD, const double omega, const double lambdas);

double ddpwfdptsl_dCbar(const void *parameters, double u);

double dpwfdptsl(const modelparameters *p, double t);

double ddpwfdptsl_dC(const modelparameters *p, double t);

#endif
