#ifndef DPWFFCF_H
#define DPWFFCF_H

#include "modelparam.h"

double dpwffcfbar(const void *parameters, double u);

double ddpwffcf_dCbar(const void *parameters, double u);

double dpwffcf(const modelparameters *p, double t);

double ddpwffcf_dC(const modelparameters *p, double t);

#endif
