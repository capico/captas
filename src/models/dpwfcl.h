#ifndef DPWFCL_H
#define DPWFCL_H

#include "modelparam.h"

double dpwfcl(const modelparameters *p, double t);

double ddpwfcl_dk(const modelparameters *p, double t);

double ddpwfcl_dS(const modelparameters *p, double t);

double ddpwfcl_dw1(const modelparameters *p, double t);

double ddpwfcl_dw2(const modelparameters *p, double t);

#endif
