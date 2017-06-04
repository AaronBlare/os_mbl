#include "TranspFs.h"

void transpFs(Model *m)
{
  int N_mat = m->N_mat;
  crsMatrix **f_mat = m->f_mat;
  crsMatrix **f_H_mat = new crsMatrix*[N_mat];
  for(int i = 0; i < N_mat; i++)
  {
    f_H_mat[i] = new crsMatrix(*(f_mat[i]));
    Transpose(*(f_mat[i]),*(f_H_mat[i]), false);
//    printMatrixVal(f_H_mat[i]);
  }
  m->f_H_mat = f_H_mat;
}

