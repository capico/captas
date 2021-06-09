#ifndef MODELS_H
#define MODELS_H

/*
Models:
dpwf        for homogeneous infinite (S and C)
dpwft       for homogeneous infinite, transformed parameters (S and C)
dpwfnfob    for homogeneous circular no-flow outer boundary (S and C)
dpwfcpob    for homogeneous circular constant pressure outer boundary (S and C)
dpwff       for homogeneous reservoir, linear sealing fault (S and C)
dpwfc       for homogeneous reservoir, parallel linear sealing faults (S and C)
dpwfcl      for homogeneous reservoir, parallel linear sealing (S and C = 0)
dpwfrect    for homogeneous rectangular reservoir (S, C)
*/

#include "dpwf.h"
#include "dpwft.h"
#include "dpwfnfob.h"
#include "dpwfcpob.h"
#include "dpwff.h"
#include "dpwfc.h"
#include "dpwfcl.h"
#include "dpwfrect.h"

#endif