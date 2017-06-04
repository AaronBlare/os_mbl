#include "linSolv.h"
#include "calcODE.h"
#include <math.h>
#include "mkl.h"

void linSolv(Model *m)
{
  crsMatrix * Gs = m->Gs;
  dcomplex  * Ks = m->Ks;
  int N_mat = m->N_mat;
  dcomplex  * RhoF = m->RhoF;
  
  for (int i = 0; i < Gs->NZ; i++)
    Gs->Col[i]++;
  for (int j = 0; j <= N_mat; j++)
  {
    Gs->RowIndex[j]++;
  }


  MKL_INT mtype = 13; /* Real unsymmetric matrix */
  
  MKL_INT nrhs = 1; /* Number of right hand sides. */
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  void *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  MKL_INT i;
  double ddum; /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  iparm[0] = 1; /* No solver default */
  iparm[1] = 2; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm[2] = 1;
  iparm[3] = 0; /* No iterative-direct algorithm */
  iparm[4] = 0; /* No user fill-in reducing permutation */
  iparm[5] = 0; /* Write solution into x */
  iparm[6] = 0; /* Not in use */
  iparm[7] = 2; /* Max numbers of iterative refinement steps */
  iparm[8] = 0; /* Not in use */
  iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0; /* Not in use */
  iparm[12] = 0; /* Not in use */
  iparm[13] = 0; /* Output: Number of perturbed pivots */
  iparm[14] = 0; /* Not in use */
  iparm[15] = 0; /* Not in use */
  iparm[16] = 0; /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0; /* Output: Numbers of CG Iterations */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 0; /* Print statistical information in file */
  error = 0; /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    pt[i] = 0;
  }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, Gs->Value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    exit(1);
  }
//  printf("\nReordering completed ... ");
//  printf("\nNumber of nonzeros in factors = %d", iparm[17]);
//  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, Gs->Value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    exit(2);
  }
//  printf("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
  phase = 33;
  iparm[7] = 2; /* Max numbers of iterative refinement steps. */
  /* Set right hand side to one. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, Gs->Value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, Ks, RhoF, &error);
  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    exit(3);
  }
  
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
  phase = -1; /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, &ddum, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);

  for (int i = 0; i < Gs->NZ; i++)
    Gs->Col[i]--;
  for (int j = 0; j <= N_mat; j++)
  {
    Gs->RowIndex[j]--;
  }
}


void linSolvCheck(Model *m)
{
  crsMatrix * Gs = m->Gs;
  dcomplex  * Ks = m->Ks;
  int N_mat = m->N_mat;
  dcomplex  * RhoF = m->RhoF;
  
  dcomplex * Ks_r = new dcomplex[N_mat];

  multMatVec(Gs, RhoF, Ks_r);

  double m1, m2;

  m1 = fabs(Ks[0].re - Ks_r[0].re);
  m2 = fabs(Ks[0].im - Ks_r[0].im);

  for(int i = 1; i < N_mat; i++)
  {
    if(m1 < fabs(Ks[0].re - Ks_r[0].re)) m1 = fabs(Ks[0].re - Ks_r[0].re);
    if(m2 < fabs(Ks[0].im - Ks_r[0].im)) m2 = fabs(Ks[0].im - Ks_r[0].im);
  }

  printf("diff (%1.16lf + i %1.16lf) \n", m1, m2);

  delete [] Ks_r;
}

void linSolvReal(Model *m)
{
  crsMatrix * Gs = m->Gs;
  dcomplex  * Ks = m->Ks;
  int N_mat = m->N_mat;
  dcomplex  * RhoF = m->RhoF;
  double * value, *res;
  value = new double[Gs->NZ];
  res = new double[Gs->N];

  for (int i = 0; i < Gs->NZ; i++)
  {
    Gs->Col[i]++;
    value[i] = Gs->Value[i].re;
  }
  for (int j = 0; j <= N_mat; j++)
  {
    Gs->RowIndex[j]++;
  }


  MKL_INT mtype = 11; /* Real unsymmetric matrix */
  
  MKL_INT nrhs = 1; /* Number of right hand sides. */
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  void *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  MKL_INT i;
  double ddum; /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  iparm[0] = 1; /* No solver default */
  iparm[1] = 2; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm[2] = 1;
  iparm[3] = 0; /* No iterative-direct algorithm */
  iparm[4] = 0; /* No user fill-in reducing permutation */
  iparm[5] = 0; /* Write solution into x */
  iparm[6] = 0; /* Not in use */
  iparm[7] = 2; /* Max numbers of iterative refinement steps */
  iparm[8] = 0; /* Not in use */
  iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0; /* Not in use */
  iparm[12] = 0; /* Not in use */
  iparm[13] = 0; /* Output: Number of perturbed pivots */
  iparm[14] = 0; /* Not in use */
  iparm[15] = 0; /* Not in use */
  iparm[16] = 0; /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0; /* Output: Numbers of CG Iterations */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 0; /* Print statistical information in file */
  error = 0; /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    pt[i] = 0;
  }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    exit(1);
  }
//  printf("\nReordering completed ... ");
//  printf("\nNumber of nonzeros in factors = %d", iparm[17]);
//  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    exit(2);
  }
//  printf("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
  phase = 33;
  iparm[7] = 2; /* Max numbers of iterative refinement steps. */
  /* Set right hand side to one. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, Ks, res, &error);
  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    exit(3);
  }
  
  for(int i = 0; i < Gs->N; i++)
  {
    RhoF[i].re = res[i];
    RhoF[i].im = 0.0;
  }

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
  phase = -1; /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &N_mat, &ddum, Gs->RowIndex, Gs->Col, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);

  for (int i = 0; i < Gs->NZ; i++)
    Gs->Col[i]--;
  for (int j = 0; j <= N_mat; j++)
  {
    Gs->RowIndex[j]--;
  }
  delete []res;
  delete []value;
}