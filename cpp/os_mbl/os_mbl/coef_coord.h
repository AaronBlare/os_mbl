#ifndef __COEF_COORD__
#define __COEF_COORD__

#include "mkl.h"

struct Tensor_Coordinates
{
  MKL_Complex16 * data;
  unsigned int * coord1;
  unsigned int * coord2;
  unsigned int * coord3;
  unsigned int * hash;
  unsigned int k;
};

struct Tensor_Coordinates_1
{
  MKL_Complex16 * data;
  unsigned int coord1;
  unsigned int * coord2;
  unsigned int * coord3;
  unsigned int N;
};

void sort_matrix(Tensor_Coordinates * matrix);
void fijk_coord(Tensor_Coordinates * f_ijk, int N);
void dijk_coord(Tensor_Coordinates * d_ijk, int N);
void free_matrix(Tensor_Coordinates * matrix);

void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat);


#endif