#ifndef DPWFNFOB_H
#define DPWFNFOB_H

#include "modelparam.h"

double dpwfnfobbar(const void *parameters, double u);

double ddpwfnfob_dCbar(const void *parameters, double u);

double ddpwfnfob_dSbar(const void *parameters, double u);

double ddpwfnfob_drebar(const void *parameters, double u);

double ddpwfnfob_dkbar(const void *parameters, double u);

double dpwfnfob(const modelparameters *p, double t);

double ddpwfnfob_dC(const modelparameters *p, double t);

double ddpwfnfob_dS(const modelparameters *p, double t);

double ddpwfnfob_dre(const modelparameters *p, double t);

double ddpwfnfob_dk(const modelparameters *p, double t);

#endif
