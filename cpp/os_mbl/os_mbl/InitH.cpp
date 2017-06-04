#include "initH.h"
#include <math.h>
#include <mkl.h>
#include "data.h"


crsMatrix * create_HE_matrix(Model * m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * HE = new crsMatrix(N + 1, N + 1);
	int * RowIndex = HE->RowIndex;
	int * Col = HE->Col;
	dcomplex * Value = HE->Value;

	double * diag_array = new double [N + 1];

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		vector<int> vb = convert_int_to_vector_of_bits(cd.id_to_x[state_id], cd.Nc);
		double sum = 0.0;
		for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
		{
			sum += double(vb[cell_id]) * cd.energies[cell_id];
		}
		sum *= 2.0 * cp.W;

		diag_array[state_id] = sum;
	}

	int i, j;
	RowIndex[0] = 0; 
	for(i = 0; i < N+1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for(i = 0; i < N + 1; i++)
	{
		Value[i].re = diag_array[i];
	}

	delete[] diag_array;

	return HE;
}

void init_hE_vector(Model * m, ConfigData &cd, ConfigParam &cp)
{
	crsMatrix * HE = create_HE_matrix(m, cd, cp), *res;
	m->HE = HE;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	dcomplex * hE = m->hE;

	int i = 0;
	for(i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*HE, *(Fs->F[i+1]), *res);
		hE[i] = trace(*res);

		delete res;
	}
}

crsMatrix * create_HU_matrix(Model * m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * HU = new crsMatrix(N + 1, N + 1);
	int * RowIndex = HU->RowIndex;
	int * Col = HU->Col;
	dcomplex * Value = HU->Value;

	double * hu_diag = new double [N + 1];

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		hu_diag[state_id] = cp.U * bit_count(cd.id_to_x[state_id] & (cd.id_to_x[state_id] << 1));
	}

	int i, j;
	RowIndex[0] = 0; 
	for(i = 0; i < N+1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for(i = 0; i < N + 1; i++)
	{
		Value[i].re = hu_diag[i];
	}

	delete [] hu_diag;

	return HU;
}

void init_hU_vector(Model * m, ConfigData &cd, ConfigParam &cp)
{
	crsMatrix * HU = create_HU_matrix(m, cd, cp), *res;
	m->HU = HU;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	dcomplex * hU = m->hU;

	int i = 0;
	for(i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*HU, *(Fs->F[i+1]), *res);
		hU[i] = trace(*res);

		delete res;
	}
}

crsMatrix * create_HJ_matrix(Model * m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * HJ;

	if (m->conf.bc == 1)
	{
		vector<int> cols;
		int * num_in_rows = new int[cd.Ns];
		int num_elems = 0;
		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			num_in_rows[state_id_1] = 0;
		}

		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
			{
				if (cd.adjacement[cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2]] > 0)
				{
					num_elems++;
					num_in_rows[state_id_1]++;
					cols.push_back(state_id_2);
				}
			}
		}

		HJ = new crsMatrix(N + 1, num_elems);
		int * RowIndex = HJ->RowIndex;
		int * Col = HJ->Col;
		dcomplex * Value = HJ->Value;

		RowIndex[0] = 0;
		for (int row_id = 1; row_id <= cd.Ns; row_id++)
		{
			RowIndex[row_id] = RowIndex[row_id-1] + num_in_rows[row_id-1];
		}

		for (int nz_id = 0; nz_id < num_elems; nz_id++)
		{
			Col[nz_id] = cols[nz_id];
			Value[nz_id].re = -cp.J;
		}

		delete[] num_in_rows;
		cols.clear();
	}
	else
	{
		vector<int> cols;
		int * num_in_rows = new int[cd.Ns];
		int num_elems = 0;
		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			num_in_rows[state_id_1] = 0;
		}

		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
			{
				if (cd.adjacement[cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2]] > 0)
				{
					num_elems++;
					num_in_rows[state_id_1]++;
					cols.push_back(state_id_2);
				}
			}
		}

		HJ = new crsMatrix(N + 1, num_elems);
		int * RowIndex = HJ->RowIndex;
		int * Col = HJ->Col;
		dcomplex * Value = HJ->Value;

		RowIndex[0] = 0;
		for (int row_id = 1; row_id <= cd.Ns; row_id++)
		{
			RowIndex[row_id] = RowIndex[row_id - 1] + num_in_rows[row_id - 1];
		}

		for (int nz_id = 0; nz_id < num_elems; nz_id++)
		{
			Col[nz_id] = cols[nz_id];
			Value[nz_id].re = -cp.J;
		}

		delete[] num_in_rows;
		cols.clear();
	}

	return HJ;
}

void init_hJ_vector(Model * m, ConfigData &cd, ConfigParam &cp)
{
	crsMatrix * HJ = create_HJ_matrix(m, cd, cp), *res;
	m->HJ = HJ;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	dcomplex * hJ = m->hJ;

	int i = 0;
	for(i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*HJ, *(Fs->F[i+1]), *res);
		hJ[i] = trace(*res);

		delete res;
	}
}

void init_H0(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * H0 = new crsMatrix();

	crsMatrix * subSum1 = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->HE), sum, *(m->HU), *subSum1);
	SparseMKLAdd(*subSum1, sum, *(m->HJ), *H0);

	m->H0 = H0;
}
