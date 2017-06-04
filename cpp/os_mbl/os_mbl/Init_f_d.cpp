#include "Init_f_d.h"
#include "coef_coord.h"

#include <vector>

using namespace std;

void init_f_d(Model *m)
{
  int N = m->N;
  int N_mat = m->N_mat;

  FMatrixs * Fs = m->Fs;

  crsMatrix * Fjk, *Fkj, *Fk, *Fsum, *Fsub, *Fres;

  dcomplex add, sub;
  add.re =  1.0; add.im = 0.0;
  sub.re = -1.0; sub.im = 0.0;

  vector<pair<int , dcomplex> > * f_std;
  vector<pair<int , dcomplex> > * d_std;
  
  m->d_mat = new crsMatrix * [N_mat];
  m->f_mat = new crsMatrix * [N_mat];

  int i, j, k, cnt;

  Fjk = new crsMatrix;
  Fkj = new crsMatrix;

  Fsum = new crsMatrix;
  Fsub = new crsMatrix;
        
  Fres = new crsMatrix;

  for( i = 0; i < N_mat; i++)
  {
    f_std = new vector<pair<int , dcomplex> >[N_mat];
    d_std = new vector<pair<int , dcomplex> >[N_mat];

    for( j = 0; j < N_mat; j++)
    {
      for( k = 0; k < N_mat; k++)
      {
        
        Fk  = new crsMatrix(*(Fs->F[k + 1]));
        
        SparseMKLMult(*(Fs->F[j + 1]), *Fk, *Fjk, true);
        SparseMKLMult(*Fk, *(Fs->F[j + 1]), *Fkj, true);
        delete Fk;

        SparseMKLAdd(*Fjk, add, *Fkj, *Fsum, true);
        SparseMKLAdd(*Fjk, sub, *Fkj, *Fsub, true);

        
        SparseMKLMult(*(Fs->F[i + 1]), *Fsub, *Fres, true);
        cnt = trace_struct(*Fres);
        if(cnt > 0)
        {
          dcomplex val, tr;
          tr = trace(*Fres);
          val.re = tr.im;
          val.im = -tr.re;
          f_std[j].push_back( make_pair(k, val));
        }
        
        SparseMKLMult(*(Fs->F[i + 1]), *Fsum, *Fres, true);
        cnt = trace_struct(*Fres);
        if(cnt > 0)
        {
          dcomplex val;
          val = trace(*Fres);
          d_std[j].push_back( make_pair(k, val));
        }
      }
    }

    m->f_mat[i] = stdToCrs(f_std, N_mat);
    m->d_mat[i] = stdToCrs(d_std, N_mat);

//    printMatrixVal(m->f_mat[i]);
//    printMatrixVal(m->d_mat[i]);

    delete [] f_std;
    delete [] d_std;
  }
  delete Fjk;
  delete Fkj;

  delete Fsum;
  delete Fsub;

  delete Fres;

}

crsMatrix * TensorToCrs(Tensor_Coordinates &t_ijk, int s, int f, int N)
{
   crsMatrix *  mat = new crsMatrix(N, f - s);
   int i = 0;
   dcomplex val;
   
   for(i = 0; i < f - s; i++)
   {
     val.re = t_ijk.data[s + i].real;
     val.im = t_ijk.data[s + i].imag;
     mat->Value[i] = val;
     mat->Col  [i] = t_ijk.coord3[s + i];
   }

   for(i = 0; i < N + 1; i++)
   {
     mat->RowIndex[i] = 0;
   }

   for(i = 0; i < f - s; i++)
   {
     mat->RowIndex[t_ijk.coord2[s + i] + 1]++;
   }

   int cnt = 0;

   for(i = 0; i < N; i++)
   {
     mat->RowIndex[i] = cnt;
     cnt += mat->RowIndex[i + 1];
   }
   mat->RowIndex[i] = cnt;

   return mat;
}

crsMatrix * TensorToCrs(Tensor_Coordinates_1 &t_ijk, int s, int f, int N)
{
   crsMatrix *  mat = new crsMatrix(N, f - s);
   int i = 0;
   dcomplex val;
   
   for(i = 0; i < f - s; i++)
   {
     val.re = t_ijk.data[s + i].real;
     val.im = t_ijk.data[s + i].imag;
     mat->Value[i] = val;
     mat->Col  [i] = t_ijk.coord3[s + i];
   }

   for(i = 0; i < N + 1; i++)
   {
     mat->RowIndex[i] = 0;
   }

   for(i = 0; i < f - s; i++)
   {
     mat->RowIndex[t_ijk.coord2[s + i] + 1]++;
   }

   int cnt = 0;

   for(i = 0; i < N; i++)
   {
     mat->RowIndex[i] = cnt;
     cnt += mat->RowIndex[i + 1];
   }
   mat->RowIndex[i] = cnt;

   return mat;
}

void init_f_d_valentin(Model *m)
{
  int N = m->N;
  int N_mat = m->N_mat;
  int i;

  m->d_mat = new crsMatrix * [N_mat];
  m->f_mat = new crsMatrix * [N_mat];

//  vector<pair<int , dcomplex> > * f_std;
//  vector<pair<int , dcomplex> > * d_std;

  Tensor_Coordinates f_ijk;
  Tensor_Coordinates d_ijk;
  Tensor_Coordinates_1 * f_1 = new Tensor_Coordinates_1[N_mat];
  Tensor_Coordinates_1 * d_1 = new Tensor_Coordinates_1[N_mat];
  
  fijk_coord(&f_ijk, N+1);
  sort_matrix(&f_ijk, f_1, N_mat);
//  sort_matrix(&f_ijk);

  dijk_coord(&d_ijk, N+1);
  sort_matrix(&d_ijk, d_1, N_mat);
//  sort_matrix(&d_ijk);

  int k, s, f;

//  k = 0;
  for( i = 0; i < N_mat; i++)
  {
//    s = k;
//    while(f_ijk.coord1[k] == i)
//    {
//      k++;
//    }
//    f = k;
//    m->f_mat[i] = TensorToCrs(f_ijk, s, f, N_mat);
    m->f_mat[i] = TensorToCrs(f_1[i], 0, f_1[i].N, N_mat);
  }

//  k = 0;
  for( i = 0; i < N_mat; i++)
  {
//    s = k;
//    while(d_ijk.coord1[k] == i)
//    {
//      k++;
//    }
//    f = k;
//    m->d_mat[i] = TensorToCrs(d_ijk, s, f, N_mat);
    m->d_mat[i] = TensorToCrs(d_1[i], 0, d_1[i].N, N_mat);
  }

  for( i = 0; i < N_mat; i++)
  {
    delete[] f_1[i].coord2;
    delete[] f_1[i].coord3;
    delete[] f_1[i].data;
    delete[] d_1[i].coord2;
    delete[] d_1[i].coord3;
    delete[] d_1[i].data;
  }
  delete[] f_1;
  delete[] d_1;

  free_matrix(&f_ijk);
  free_matrix(&d_ijk);
}