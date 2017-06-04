#include "CalcRho.h"

void calcRho(Model *m)
{
  FMatrixs  * Fs = m->Fs;
  crsMatrix * Rho = new crsMatrix(*(Fs->F[0]));
  crsMatrix * Rho_tmp = new crsMatrix(*(Fs->F[0]));
  crsMatrix * tmp;
  dcomplex  * RhoF = m->RhoF;
  int N_mat = m->N_mat;
  int N = m->N;
 
  int i;
  for(i = 0; i < Rho->NZ; i++)
  {
    Rho->Value[i].re /= N + 1;
    Rho->Value[i].im /= N + 1;
  }

  for(int i = 0; i < N_mat; i++)
  {
    SparseMKLAdd(*Rho, RhoF[i], *(Fs->F[i + 1]), *Rho_tmp, true); 
    tmp = Rho; Rho = Rho_tmp; Rho_tmp = tmp;
  }
  delete Rho_tmp;

  m->Rho = Rho;
}

