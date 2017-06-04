#include "CalcTraseRO2.h"
#include "Matrix_op.h"
#include <stdio.h>

void calcTraseRO2(Model *m)
{
  crsMatrix * Rho = m->Rho;
  crsMatrix * Rho2 = new crsMatrix;
  
  int i, j;
  toOneBase(*Rho);

  SparseMKLMultOne(*Rho, *Rho, *Rho2);

  toZeroBase(*Rho);
  toZeroBase(*Rho2);
  
  dcomplex tr= trace(*Rho2);
  
  delete Rho2;
  
  char fileName[512];
  sprintf(fileName, "traseRho2.txt");
  FILE * f = fopen(fileName, "w");
  
  fprintf(f, "%0.16lf %0.16lf\n", tr.re, tr.im);
  
  fclose(f);
}
