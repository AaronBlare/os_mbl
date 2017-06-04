#include "coef_coord.h"
#include "stdio.h"
#include "math.h"
#include "string.h"


#define IND(i, j, k) (i)*(N * N - 1)*(N * N - 1) + (j)*(N * N - 1) + (k)

#define IndS(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2)
#define IndJ(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2 + 1)
#define IndD(l)    (N * (N-1) + l - 1)


void sort_matrix(Tensor_Coordinates * matrix)
{
  unsigned int key;
  unsigned int pos;
  MKL_Complex16 data_tmp;
  unsigned int int_tmp;
  
  for(unsigned int i = 0; i < matrix[0].k; i++)
  {
    key = matrix[0].hash[i];
    pos = i;
    for(unsigned int j = i + 1; j < matrix[0].k; j++)
    {
      if (matrix[0].hash[j] < key)
      {
        key = matrix[0].hash[j];
        pos = j;
      }
    }
    if(pos != i)
    {
      data_tmp = matrix[0].data[i];
      matrix[0].data[i] = matrix[0].data[pos];
      matrix[0].data[pos] = data_tmp;

      int_tmp = matrix[0].coord1[i];
      matrix[0].coord1[i] = matrix[0].coord1[pos];
      matrix[0].coord1[pos] = int_tmp;

      int_tmp = matrix[0].coord2[i];
      matrix[0].coord2[i] = matrix[0].coord2[pos];
      matrix[0].coord2[pos] = int_tmp;

      int_tmp = matrix[0].coord3[i];
      matrix[0].coord3[i] = matrix[0].coord3[pos];
      matrix[0].coord3[pos] = int_tmp;

      int_tmp = matrix[0].hash[i];
      matrix[0].hash[i] = matrix[0].hash[pos];
      matrix[0].hash[pos] = int_tmp;
    }
  }
}

void allocMemMat(Tensor_Coordinates_1 * mat)
{
  mat->coord2 = new unsigned  int[mat->N];
  mat->coord3 = new unsigned  int[mat->N];
  mat->data   = new MKL_Complex16[mat->N];
}

void swap_row(Tensor_Coordinates_1 &mat, int i, int j)
{
  unsigned int t;
  MKL_Complex16 v;
  t = mat.coord2[i]; mat.coord2[i] = mat.coord2[j]; mat.coord2[j] = t; 
  t = mat.coord3[i]; mat.coord3[i] = mat.coord3[j]; mat.coord3[j] = t;
  v = mat.data  [i]; mat.data  [i] = mat.data  [j]; mat.data  [j] = v;
}

void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat)
{
  int i = 0;
  int NN = matrix->k;
  for(i = 0; i < Nmat; i++)
  {
    mat_res[i].coord1 = i;
    mat_res[i].N = 0;
  }

  for(i = 0; i < NN; i++)
  {
    int jj = matrix->coord1[i];
    if(jj < Nmat)
    {
      mat_res[jj].N++;
    }
    else 
    {
      printf("*");
    }
  }

  for(i = 0; i < Nmat; i++)
  {
    allocMemMat(mat_res+i);
    mat_res[i].N = 0;
  }

  int jj, ii, k;
  for(i = 0; i < NN; i++)
  {
    jj = matrix->coord1[i];
    ii = mat_res[jj].N;
    mat_res[jj].coord2[ii] = matrix->coord2[i];
    mat_res[jj].coord3[ii] = matrix->coord3[i];
    mat_res[jj].data  [ii] = matrix->data  [i];
    mat_res[jj].N++;
  }

  for(i = 0; i < Nmat; i++)
  {
    Tensor_Coordinates_1 tmp = mat_res[i];
    for(ii = 0; ii < tmp.N - 1; ii++)
    {
      for(jj = 0; jj < tmp.N - 1; jj++)
      {
        if(tmp.coord2[jj] > tmp.coord2[jj+1])
        {
          swap_row(tmp, jj, jj + 1);
        }
        else
        {
          if(tmp.coord2[jj] == tmp.coord2[jj+1])
          {
            if(tmp.coord3[jj] > tmp.coord3[jj+1])
            {
              swap_row(tmp, jj, jj + 1);
            }
          }
        }
      }
    }
  }
}

void fijk_coord(Tensor_Coordinates * f_ijk, int N)
{
  unsigned int size = 5 * N * N * N - 9 * N * N - 2 * N + 6;
  f_ijk[0].data = new MKL_Complex16[size];
  memset(f_ijk[0].data, 0, size * sizeof(MKL_Complex16));
  f_ijk[0].coord1 = new unsigned int[size];
  memset(f_ijk[0].coord1, 0, size * sizeof(unsigned int));
  f_ijk[0].coord2 = new unsigned int[size];
  memset(f_ijk[0].coord2, 0, size * sizeof(unsigned int));
  f_ijk[0].coord3 = new unsigned int[size];
  memset(f_ijk[0].coord3, 0, size * sizeof(unsigned int));
  f_ijk[0].hash = new unsigned int[size];
  memset(f_ijk[0].hash, 0, size * sizeof(unsigned int));
  f_ijk[0].k = 0;

  for(int i = 1; i < N; i++)  
  {
    for(int j = i + 1; j < N; j++)
    {
      f_ijk[0].data[f_ijk[0].k].real = (i)/sqrt((double)(i) * (i + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndD(i); 
      f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndD(i), IndJ(i, j), IndS(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndD(i); 
      f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndD(i), IndS(i, j), IndJ(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndD(i);
      f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(i), IndS(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = (i)/sqrt((double)(i) * (i + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndD(i);
      f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(i), IndJ(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = (i)/sqrt((double)(i) * (i + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndD(i);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(i));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndD(i);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(i));
      f_ijk[0].k++;
    }
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = i + 1; j < N; j++)
    {
      for(int m = i + 1; m < j; m++)
      {
        f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt((double)(m) * (m + 1));
        f_ijk[0].coord1[f_ijk[0].k] = IndD(m); 
        f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
        f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
        f_ijk[0].hash[f_ijk[0].k] = IND(IndD(m), IndJ(i, j), IndS(i, j));
        f_ijk[0].k++;

        f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        f_ijk[0].coord1[f_ijk[0].k] = IndD(m); 
        f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
        f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
        f_ijk[0].hash[f_ijk[0].k] = IND(IndD(m), IndS(i, j), IndJ(i, j));
        f_ijk[0].k++;

        f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j); 
        f_ijk[0].coord2[f_ijk[0].k] = IndD(m);
        f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
        f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(m), IndS(i, j));
        f_ijk[0].k++;

        f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt((double)(m) * (m + 1));
        f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j); 
        f_ijk[0].coord2[f_ijk[0].k] = IndD(m);
        f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
        f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(m), IndJ(i, j));
        f_ijk[0].k++;

        f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt((double)(m) * (m + 1));
        f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j); 
        f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
        f_ijk[0].coord3[f_ijk[0].k] = IndD(m);
        f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(m));
        f_ijk[0].k++;

        f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j); 
        f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
        f_ijk[0].coord3[f_ijk[0].k] = IndD(m);
        f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(m));
        f_ijk[0].k++;
      }
    }
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = i + 1; j < N; j++)
    {
      f_ijk[0].data[f_ijk[0].k].real = -(1 + j)/sqrt((double)(j) * (j + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndD(j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndS(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = (1 + j)/sqrt((double)(j) * (j + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndD(j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndD(j), IndS(i, j), IndJ(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = (1 + j)/sqrt((double)(j) * (j + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndD(j);
      f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndS(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = -(1 + j)/sqrt((double)(j) * (j + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndD(j);
      f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(j), IndJ(i, j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = -(1 + j)/sqrt((double)(j) * (j + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndD(j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(j));
      f_ijk[0].k++;

      f_ijk[0].data[f_ijk[0].k].real = (1 + j)/sqrt((double)(j) * (j + 1));
      f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j); 
      f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
      f_ijk[0].coord3[f_ijk[0].k] = IndD(j);
      f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(j));
      f_ijk[0].k++;
    }
  }

  for(int i = 0; i < N-1; i++)   
  {
    for(int j = i + 1; j < N; j++)
    {
      for(int k = 0; k < N; k++)
      {
        if(k > j)
        {
          if(k > i)
          {
            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0); //-i*Jjk*[Sij,Sik]  //and symmetric Ski, Jkj
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(j, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(j, k), IndS(i, j), IndS(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0); //-i*Jik*[Sij,Sjk]  //and symmetric Skj, Jki
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(j, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndS(j, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0); //-i*Sjk*[Sij,Jik]  //and symmetric Skj, Jki
            f_ijk[0].coord1[f_ijk[0].k] = IndS(j, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(j, k), IndS(i, j), IndJ(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0); //-i*Sik*[Sij,Jjk]  //and symmetric Ski, Jkj
            f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(j, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndJ(j, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0); //-i*Sjk*[Jij,Sik]  //and symmetric Ski, Skj
            f_ijk[0].coord1[f_ijk[0].k] = IndS(j, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(j, k), IndJ(i, j), IndS(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0); //-i*Sik*[Jij,Sjk]  //and symmetric Skj, Ski
            f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(j, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndS(j, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0); //-i*Jjk*[Jij,Jik]  //and symmetric Jkj, Jki
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(j, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(j, k), IndJ(i, j), IndJ(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0); //-i*Jik*[Jij,Jjk]  //and symmetric Jki, Jkj
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(j, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndJ(j, k));
            f_ijk[0].k++;
          }
        }
        if(k < j)
        {
          if(k > i)
          {
            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndS(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndS(k, j));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndJ(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndJ(k, j));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndS(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndS(k, j));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndJ(i, k));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndJ(k, j));
            f_ijk[0].k++;
          }
          if(k < i)
          {
            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(k, i);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndS(k, i));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, i); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, i), IndS(i, j), IndS(k, j));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, i);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndJ(k, i));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(k, i); 
            f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, i), IndS(i, j), IndJ(k, j));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(k, i);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndS(k, i));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndS(k, i); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, i), IndJ(i, j), IndS(k, j));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = 1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, i);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndJ(k, i));
            f_ijk[0].k++;

            f_ijk[0].data[f_ijk[0].k].real = -1.0/sqrt(2.0);
            f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, i); 
            f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
            f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
            f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, i), IndJ(i, j), IndJ(k, j));
            f_ijk[0].k++;
          }
        }
      }
    }
  }
}

void dijk_coord(Tensor_Coordinates * d_ijk, int N)
{
  unsigned int size = 6 * N * N * N - (N * (21 * N + 7))/2 + 1;
  d_ijk[0].data = new MKL_Complex16[size];
  memset(d_ijk[0].data, 0, size * sizeof(MKL_Complex16));
  d_ijk[0].coord1 = new unsigned int[size];
  memset(d_ijk[0].coord1, 0, size * sizeof(unsigned int));
  d_ijk[0].coord2 = new unsigned int[size];
  memset(d_ijk[0].coord2, 0, size * sizeof(unsigned int));
  d_ijk[0].coord3 = new unsigned int[size];
  memset(d_ijk[0].coord3, 0, size * sizeof(unsigned int));
  d_ijk[0].hash = new unsigned int[size];
  memset(d_ijk[0].hash, 0, size * sizeof(unsigned int));
  d_ijk[0].k = 0;

  for(int i = 1; i < N; i++)  
  {
    for(int j = i + 1; j < N; j++)
    {
      d_ijk[0].data[d_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndD(i); 
      d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndS(i, j), IndS(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
      d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(i), IndS(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
      d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndJ(i, j), IndJ(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
      d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(i), IndJ(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(i));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = -(i)/sqrt((double)(i) * (i + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(i));
      d_ijk[0].k++;
    }
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = i + 1; j < N; j++)
    {
      for(int m = i + 1; m < j; m++)
      {
        d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndD(m); 
        d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndD(m), IndS(i, j), IndS(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndD(m);
        d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(m), IndS(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndD(m);
        d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndD(m), IndJ(i, j), IndJ(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndD(m);
        d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(m), IndJ(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndD(m);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(m));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt((double)(m) * (m + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndD(m);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(m));
        d_ijk[0].k++;
      }
    }
  }

  for(int j = 2; j < N; j++)
  {
    int i = 0;
    d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndD(j); 
    d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
    d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndS(i, j), IndS(i, j));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
    d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
    d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(j), IndS(i, j));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
    d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
    d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
    d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
    d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
    d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
    d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(j));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
    d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
    d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
    d_ijk[0].k++;
  }

  for(int i = 1; i < N; i++)
  {
    for(int j = i + 1; j < N; j++)
    {
      d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndD(j); 
      d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndS(i, j), IndS(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
      d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(j), IndS(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
      d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
      d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = (1 - j)/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
      d_ijk[0].k++;
    }
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = i + 1; j < N; j++)
    {
      for(int z = j + 1; z < N; z++)
      {
        d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(z) * (z + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndD(z); 
        d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndD(z), IndS(i, j), IndS(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(z) * (z + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndD(z);
        d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(z), IndS(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(z) * (z + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndD(z);
        d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndD(z), IndJ(i, j), IndJ(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(z) * (z + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndD(z);
        d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(z), IndJ(i, j));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(z) * (z + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndD(z);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(z));
        d_ijk[0].k++;

        d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(z) * (z + 1));
        d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
        d_ijk[0].coord3[d_ijk[0].k] = IndD(z);
        d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(z));
        d_ijk[0].k++;
      }
    }
  }

  for(int i = 0; i < N-1; i++)   
  {
    for(int j = i + 1; j < N; j++)
    {
      for(int k = 0; k < N; k++)
      {
        if(k > j)
        {
          if(k > i)
          {
            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);  //Sjk*{Sij,Sik}  //and symmetric Ski, Skj
            d_ijk[0].coord1[d_ijk[0].k] = IndS(j, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(j, k), IndS(i, j), IndS(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);  //Sik*{Sij,Sjk}  //and symmetric Skj, Ski
            d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(j, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndS(j, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);  //Jjk*{Sij,Jik}  //and symmetric Jkj, Jki
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(j, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(j, k), IndS(i, j), IndJ(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);  //Jik*{Sij,Jjk}  //and symmetric Jki, Jkj
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(j, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndJ(j, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = -1.0/sqrt(2.0);  //Jjk*{Jij,Sik}  //and symmetric Ski, Jkj
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(j, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(j, k), IndJ(i, j), IndS(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);  //Jik*{Jij,Sjk}  //and symmetric Skj, Jki
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(j, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndS(j, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);  //Sjk*{Jij,Jik}  //and symmetric Skj, Jki
            d_ijk[0].coord1[d_ijk[0].k] = IndS(j, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(j, k), IndJ(i, j), IndJ(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = -1.0/sqrt(2.0);  //Sik*{Jij,Jjk}  //and symmetric Ski, Jkj
            d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(j, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndJ(j, k));
            d_ijk[0].k++;
          }
        }
        if(k < j)
        {
          if(k > i)
          {
            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndS(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndS(k, j));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = -1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndJ(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = -1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndJ(k, j));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndS(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndS(k, j));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndJ(i, k));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndJ(k, j));
            d_ijk[0].k++;
          }
          if(k < i)
          {
            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(k, i);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndS(k, i));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(k, i); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, i), IndS(i, j), IndS(k, j));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, i);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndJ(k, i));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, i); 
            d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, i), IndS(i, j), IndJ(k, j));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(k, i);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndS(k, i));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = -1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, i); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, i), IndJ(i, j), IndS(k, j));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = -1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, i);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndJ(k, i));
            d_ijk[0].k++;

            d_ijk[0].data[d_ijk[0].k].real = 1.0/sqrt(2.0);
            d_ijk[0].coord1[d_ijk[0].k] = IndS(k, i); 
            d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
            d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
            d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, i), IndJ(i, j), IndJ(k, j));
            d_ijk[0].k++;
          }
        }
      }
    }
  }

  
  for(int j = 2; j < N; j++)
  {
    int i = 1;
    d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(j) * (j + 1));  
    d_ijk[0].coord1[d_ijk[0].k] = IndD(j); 
    d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
    d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndD(i), IndD(i));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndD(i); 
    d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
    d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(j));
    d_ijk[0].k++;

    d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(j) * (j + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndD(i); 
    d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
    d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(j), IndD(i));
    d_ijk[0].k++;
  }

  for(int i = 2; i < N; i++)
  {
    d_ijk[0].data[d_ijk[0].k].real = 2.0 * (1 - i)/sqrt((double)(i) * (i + 1));
    d_ijk[0].coord1[d_ijk[0].k] = IndD(i); 
    d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
    d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
    d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(i));
    d_ijk[0].k++;

    for(int j = i + 1; j < N; j++)
    {
      d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(j) * (j + 1));  
      d_ijk[0].coord1[d_ijk[0].k] = IndD(j); 
      d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndD(i), IndD(i));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndD(i); 
      d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(j));
      d_ijk[0].k++;

      d_ijk[0].data[d_ijk[0].k].real = 2.0/sqrt((double)(j) * (j + 1));
      d_ijk[0].coord1[d_ijk[0].k] = IndD(i); 
      d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
      d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
      d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(j), IndD(i));
      d_ijk[0].k++;
    }
  }
}

void print_matrix(Tensor_Coordinates * matrix, FILE * f)
{
  for (unsigned int i = 0; i < matrix[0].k; i++)
  {
    fprintf(f,"%2.3i %2.3i %2.3i %lf\n", matrix[0].coord1[i] + 1, matrix[0].coord2[i] + 1, matrix[0].coord3[i] + 1, matrix[0].data[i].real);
  }
}

void free_matrix(Tensor_Coordinates * matrix)
{
  delete(matrix[0].data);
  delete(matrix[0].coord1);
  delete(matrix[0].coord2);
  delete(matrix[0].coord3);
  delete(matrix[0].hash);
  matrix[0].k=0;
}
/*
void main()
{
  int N=11;

  Tensor_Coordinates f_ijk;
  Tensor_Coordinates d_ijk;

  fijk_coord(&f_ijk, N);
  sort_matrix(&f_ijk);

  dijk_coord(&d_ijk, N);
  sort_matrix(&d_ijk);

  FILE * file1 = fopen("fijk.txt", "w");

  print_matrix(&f_ijk, file1);

  FILE * file2 = fopen("dijk.txt", "w");

  print_matrix(&d_ijk, file2);
  
  fclose(file1);
  fclose(file2);
  free_matrix(&f_ijk);
  free_matrix(&d_ijk);
}*/