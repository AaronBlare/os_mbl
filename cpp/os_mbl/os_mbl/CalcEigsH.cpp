#include "CalcGs.h"
#include <string.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>

void calcEigsH0(Model *m)
{
	crsMatrix * H0 = m->H0;
	int N_mat = H0->N;
	dcomplex * fillMat = new dcomplex[N_mat * N_mat];
	dcomplex * DP = new dcomplex[N_mat];
	dcomplex * evl = new dcomplex[N_mat * N_mat];
	dcomplex * evr = new dcomplex[N_mat * N_mat];
	int i, j, k, s, f;

	memset(fillMat, 0, sizeof(dcomplex) * N_mat * N_mat);
	for(i = 0; i < N_mat; i++)
	{
		s = H0->RowIndex[i];
		f = H0->RowIndex[i + 1];
		for(k = s; k < f; k++)
		{
			j = H0->Col[k];
			fillMat[i * N_mat + j].re = H0->Value[k].re;
			fillMat[i * N_mat + j].im = H0->Value[k].im;
		}
	}

	double l_time = omp_get_wtime();
	int info;
	info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', N_mat, (MKL_Complex16 *)fillMat, N_mat, 
		(MKL_Complex16 *)DP, (MKL_Complex16 *) evl, 
		N_mat, (MKL_Complex16 *)evr , N_mat);
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	l_time = omp_get_wtime()-l_time;
	printf("LAPACKE_zgeev_H0 time:%lf\n", l_time);

	char fileName[512];

	sprintf(fileName, "evals.txt");

	FILE * file = fopen(fileName, "w");

	for(i = 0; i < N_mat; i++)
	{
		fprintf(file, "%1.16lf %1.16lf\n", DP[i].re, DP[i].im);
	}

	fclose(file);

	sprintf(fileName, "evecs_right_H0.txt");
	file = fopen(fileName, "w");

	for(int i = 0; i < N_mat * N_mat; i++)
	{

		fprintf(file, "%1.16lf %1.16lf\n", evr[i].re, evr[i].im);
	}

	fclose(file);

	delete [] DP;
	delete [] fillMat;

	delete [] evl;
	delete [] evr;
}