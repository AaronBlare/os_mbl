#ifndef __MATRIX__
#define __MATRIX__

#include <stdio.h>
// Complex (double is the base datatype for real and imaginary parts)
typedef struct dcomplex
{
	double re;
	double im;
} dcomplex;


struct crsMatrix
{
protected:
	int sizeMem;

public:
	int N;
	int NZ;

	dcomplex* Value;
	int* Col;
	int* RowIndex;

	crsMatrix()
	{
		N = 0;
		sizeMem = NZ = 0;
		Value = NULL;
		Col = NULL;
		RowIndex = NULL;
	}

	void resize(int _N, int _NZ = 0)
	{
		if (N < _N)
		{
			if (RowIndex != NULL)
			{
				delete[] RowIndex;
			}
			N = _N;
			RowIndex = new int[N + 1];
		}
		N = _N;

		NZ = _NZ;
		if (_NZ > sizeMem)
		{
			sizeMem = NZ;
			if (Col != NULL)
			{
				delete[] Col;
			}
			Col = new int[NZ];

			if (Value != NULL)
			{
				delete[] Value;
			}
			Value = new dcomplex[NZ];
		}
	}

	crsMatrix(int _N, int _NZ)
	{
		N = _N;
		sizeMem = NZ = _NZ;

		Value = new dcomplex[NZ];
		Col = new int[NZ];
		RowIndex = new int[N + 1];

		for (int i = 0; i < NZ; i++)
		{
			Value[i].re = 0.0;
			Value[i].im = 0.0;
		}
	}

	crsMatrix(const crsMatrix &mat)
	{
		N = mat.N;
		sizeMem = NZ = mat.NZ;

		Value = new dcomplex[NZ];
		Col = new int[NZ];
		RowIndex = new int[N + 1];

		for (int i = 0; i < N + 1; i++)
		{
			RowIndex[i] = mat.RowIndex[i];
		}

		for (int i = 0; i < NZ; i++)
		{
			Value[i].re = mat.Value[i].re;
			Value[i].im = mat.Value[i].im;
			Col[i] = mat.Col[i];
		}

	}

	~crsMatrix()
	{
		delete[] Value;
		delete[] Col;
		delete[] RowIndex;
	}
};

#endif