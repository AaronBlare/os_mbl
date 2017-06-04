#include "InitFs.h"
#include <stdio.h>
#include <math.h>

void outFs(FMatrixs *Fs)
{
	for (int i = 0; i < Fs->countF; i++)
	{
		printMatrixVal(Fs->F[i]);
	}
}

void initFs(FMatrixs *Fs, int N)
{
	int i, j, k;
	Fs->F[0] = createFeyeType(N + 1);
	k = 1;
	for (i = 0; i < N + 1; i++)
	{
		for (j = i + 1; j < N + 1; j++)
		{
			Fs->F[k] = createFPairTypeRe(N + 1, i, j); k++;
			Fs->F[k] = createFPairTypeIm(N + 1, i, j); k++;
		}
	}

	for (i = 0; i < N; i++)
	{
		if (k < Fs->countF)
		{
			Fs->F[k] = createLastType(N + 1, i); k++;
		}
		else
		{
			throw("error count calc (no mem)");
		}
	}

	if (k != Fs->countF)
	{
		throw("error count calc (countF > k)");
	}

	//outFs(Fs);
}

crsMatrix * createFeyeType(int N)
{
	crsMatrix * mat;
	mat = new crsMatrix(N, N);
	for (int i = 0; i < N; i++)
	{
		mat->Col[i] = i;
		mat->RowIndex[i] = i;
		mat->Value[i].re = 1.0;
	}
	mat->RowIndex[N] = N;

	return mat;
}
crsMatrix * createFPairTypeRe(int N, int i, int j)
{
	double val = 1.0 / sqrt(2.0);
	crsMatrix * mat;
	mat = new crsMatrix(N, 2);
	mat->Value[0].re = val;
	mat->Value[1].re = val;
	mat->Col[0] = j;
	mat->Col[1] = i;

	for (int ii = 0; ii < i + 1; ii++)
	{
		mat->RowIndex[ii] = 0;
	}
	for (int ii = i + 1; ii < j + 1; ii++)
	{
		mat->RowIndex[ii] = 1;
	}
	for (int ii = j + 1; ii <= N; ii++)
	{
		mat->RowIndex[ii] = 2;
	}

	return mat;
}
crsMatrix * createFPairTypeIm(int N, int i, int j)
{
	double val = -1.0 / sqrt(2.0);
	crsMatrix * mat;
	mat = new crsMatrix(N, 2);
	mat->Value[0].im = val;
	mat->Value[1].im = -val;
	mat->Col[0] = j;
	mat->Col[1] = i;

	for (int ii = 0; ii < i + 1; ii++)
	{
		mat->RowIndex[ii] = 0;
	}
	for (int ii = i + 1; ii < j + 1; ii++)
	{
		mat->RowIndex[ii] = 1;
	}
	for (int ii = j + 1; ii <= N; ii++)
	{
		mat->RowIndex[ii] = 2;
	}

	return mat;
}

crsMatrix * createLastType(int N, int i)
{
	crsMatrix * mat;
	mat = new crsMatrix(N, i + 2);
	int ii;
	double val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
	mat->RowIndex[0] = 0;
	for (ii = 0; ii <= i; ii++)
	{
		mat->Col[ii] = ii;
		mat->RowIndex[ii + 1] = mat->RowIndex[ii] + 1;
		mat->Value[ii].re = val;
	}
	mat->Col[ii] = ii;
	mat->RowIndex[ii + 1] = mat->RowIndex[ii] + 1;
	mat->Value[ii].re = -(i + 1) * val;
	ii++;

	for (; ii < N; ii++)
	{
		mat->RowIndex[ii + 1] = mat->RowIndex[ii];
	}

	return mat;
}