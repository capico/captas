#ifndef DELTAPWFF_H
#define DELTAPWFF_H

#include "modelparam.h"

double dpwffbar(const void *parameters, double u);

double ddpwff_dkbar(const void *parameters, double u);

double ddpwff_dCbar(const void *parameters, double u);

double ddpwff_dSbar(const void *parameters, double u);

double ddpwff_dLbar(const void *parameters, double u);

double dpwff(const modelparameters *p, double t);

double ddpwff_dk(const modelparameters *p, double t);

double ddpwff_dC(const modelparameters *p, double t);

double ddpwff_dS(const modelparameters *p, double t);

double ddpwff_dL(const modelparameters *p, double t);

#endif
