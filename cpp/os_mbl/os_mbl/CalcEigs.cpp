#include "CalcGs.h"
#include <string.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>

void calcEig(Model *m)
{
	crsMatrix * Gs = m->Gs;
	int N_mat = Gs->N;
	dcomplex * fillMat = new dcomplex[N_mat * N_mat];
	dcomplex * DP = new dcomplex[N_mat];
	int i, j, k, s, f;

	memset(fillMat, 0, sizeof(dcomplex) * N_mat * N_mat);
	for(i = 0; i < N_mat; i++)
	{
		s = Gs->RowIndex[i];
		f = Gs->RowIndex[i + 1];
		for(k = s; k < f; k++)
		{
			j = Gs->Col[k];
			fillMat[i * N_mat + j].re = Gs->Value[k].re;
			fillMat[i * N_mat + j].im = Gs->Value[k].im;
		}
	}

	double l_time = omp_get_wtime();
	int info;
	info = LAPACKE_zgeev( LAPACK_ROW_MAJOR, /*'V'*/'N', /*'V'*/'N', N_mat, (MKL_Complex16 *)fillMat, N_mat, 
		(MKL_Complex16 *)DP,NULL/* (MKL_Complex16 *)&LeftSigmaP[0][0]*/, 
		N_mat, NULL/*(MKL_Complex16 *)&SigmaP[0][0]*/, N_mat );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	l_time = omp_get_wtime()-l_time;
	printf("LAPACKE_zgeev time:%lf\n", l_time);

	char fileName[512];
	sprintf(fileName, "eigs.txt");
	FILE * file = fopen(fileName, "w");

	for(i = 0; i < N_mat; i++)
	{
		fprintf(file, "%1.16lf %1.16lf\n", DP[i].re, DP[i].im);
	}

	fclose(file);

	delete [] DP;
	delete [] fillMat;
}