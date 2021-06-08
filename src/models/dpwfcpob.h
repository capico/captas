#ifndef DPWFCPOB_H
#define DPWFCPOB_H

#include "modelparam.h"

double dpwfcpobbar(const void *parameters, double u);

double ddpwfcpob_dCbar(const void *parameters, double u);

double dpwfcpob(const modelparameters *p, double t);

double ddpwfcpob_dC(const modelparameters *p, double t);

#endif
