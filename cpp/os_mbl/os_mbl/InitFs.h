#ifndef __INIT_FS__
#define __INIT_FS__

#include "Model.h"

void initFs(FMatrixs *Fs, int N);
void outFs(FMatrixs *Fs);

crsMatrix * createFeyeType(int N);
crsMatrix * createFPairTypeRe(int N, int i, int j);
crsMatrix * createFPairTypeIm(int N, int i, int j);
crsMatrix * createLastType(int N, int i);
#endif