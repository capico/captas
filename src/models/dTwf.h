#ifndef DTWF_H
#define DTWF_H

#include "modelparam.h"

double dTwfbar(const void *parameters, double u);

double ddTwf_dkbar(const void *parameters, double u);

double ddTwf_dCbar(const void *parameters, double u);

double ddTwf_dSbar(const void *parameters, double u);

double dTwf(const modelparameters *T, double t);

double ddTwf_dk(const modelparameters *T, double t);

double ddTwf_dC(const modelparameters *T, double t);

double ddTwf_dS(const modelparameters *T, double t);

double dr_dpi(const modelparameters *T, double t);

#endif
