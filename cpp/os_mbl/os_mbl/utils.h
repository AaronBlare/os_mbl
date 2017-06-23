#pragma once
#include "config.h"

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

#define Error( msg ) error( msg, __FUNCTION__, __FILE__, __LINE__ )

#define PI 3.1415926535897932384626433832795

#define EPS 1.0e-12

void error(const std::string& err, const char* func, const char* file, int line);

int n_choose_k(int n, int k);

int bit_count(int value);

int bit_at(int value, int position);

void print_int_array(int * data, int N);

string file_name_suffix(ConfigParam &cp, int precision);

string file_name_suffix_zev(ConfigParam &cp, int precision);

string file_name_suffix_int(ConfigParam &cp, int precision);

string file_name_suffix_tr(ConfigParam &cp, int precision);

string file_name_suffix_qj(ConfigParam &cp, int precision);

void write_double_data(string file_name, double * data, int size, int precision, bool append);

void write_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append);

void write_sparse_complex_mtx(string file_name, crsMatrix *A, int precision, bool append);

vector<int> convert_int_to_vector_of_bits(int x, int size);

vector<int> sort_doubles_with_order(vector<double> &v);
