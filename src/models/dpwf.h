#ifndef DPWF_H
#define DPWF_H

#include "modelparam.h"

double dpwfbar(const void *parameters, double u);

double ddpwf_dkbar(const void *parameters, double u);

double ddpwf_dCbar(const void *parameters, double u);

double ddpwf_dSbar(const void *parameters, double u);

double dpwf(const modelparameters *p, double t);

double ddpwf_dk(const modelparameters *p, double t);

double ddpwf_dC(const modelparameters *p, double t);

double ddpwf_dS(const modelparameters *p, double t);

double dr_dpi(const modelparameters *p, double t);

#endif
