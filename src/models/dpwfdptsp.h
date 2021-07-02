#ifndef DPWFDPTSP_H
#define DPWFDPTSP_H

#include "modelparam.h"
#include "dpwfdp.h"

double ftsp(const double uD, const double omega, const double lambdas);

double ddpwfdptsp_dCbar(const void *parameters, double u);

double dpwfdptsp(const modelparameters *p, double t);

double ddpwfdptsp_dC(const modelparameters *p, double t);

#endif
