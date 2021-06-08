#ifndef DPWFC_H
#define DPWFC_H

#include "modelparam.h"

double dpwfcbar(const void *parameters, double u);

double ddpwfc_dkbar(const void *parameters, double u);

double ddpwfc_dCbar(const void *parameters, double u);

double ddpwfc_dSbar(const void *parameters, double u);

double ddpwfc_dw1bar(const void *parameters, double u);

double ddpwfc_dw2bar(const void *parameters, double u);

double dpwfc(const modelparameters *p, double t);

double ddpwfc_dk(const modelparameters *p, double t);

double ddpwfc_dC(const modelparameters *p, double u);

double ddpwfc_dS(const modelparameters *p, double u);

double ddpwfc_dw1(const modelparameters *p, double u);

double ddpwfc_dw2(const modelparameters *p, double u);

#endif
