#ifndef __CALC_ODE__
#define __CALC_ODE__

#include "Model.h"

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res);

void initRhoODE(Model *m);
void calcODE(Model *m, double h, int cntItr, double t = 0.0);

dcomplex calcDiffIter(Model *m);

#endif