#ifndef DPWFLE_H
#define DPWFLE_H

#include "modelparam.h"

double dpwflebar(const void *parameters, double u);

double ddpwfle_dCbar(const void *parameters, double u);

double dpwfle(const modelparameters *p, double t);

double ddpwfle_dC(const modelparameters *p, double t);

#endif
