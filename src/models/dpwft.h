#ifndef DPWFT_H
#define DPWFT_H

#include "modelparam.h"
#include "dpwf.h"

double ddpwft_dk(const modelparameters *p, double t);

double ddpwft_dC(const modelparameters *p, double t);

double ddpwft_dS(const modelparameters *p, double t);

double drt_dpi(const modelparameters *p, double t);

#endif
