#include "stdToCrs.h"

crsMatrix * stdToCrs(vector<pair<int , dcomplex> > * mat, int N)
{
  crsMatrix * res = NULL;

  int NZ = 0;
  int i, j, k, Nl;

  for(i = 0; i < N; i++)
  {
    NZ += mat[i].size();
  }
  res = new crsMatrix(N, NZ);
  k=0;
  res->RowIndex[0] = 0; 
  for(i = 0; i < N; i++)
  {
    Nl = mat[i].size();
    for(j = 0; j < Nl; j++)
    {
      res->Col[k] = mat[i][j].first;
      res->Value[k] = mat[i][j].second;
      k++;
    }
    res->RowIndex[i+1] = k;
  }

  return res;
}