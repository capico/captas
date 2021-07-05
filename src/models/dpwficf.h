#ifndef DPWFICF_H
#define DPWFICF_H

#include "modelparam.h"

double dpwficfbar(const void *parameters, double u);

double ddpwficf_dCbar(const void *parameters, double u);

double dpwficf(const modelparameters *p, double t);

double ddpwficf_dC(const modelparameters *p, double t);

#endif
