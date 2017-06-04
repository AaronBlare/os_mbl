#include "Init_a1_a2.h"
#include "config.h"
#include "data.h"
#include <math.h>

crsMatrix * create_A1_diss1_matrix(Model * m, ConfigData &cd, ConfigParam &cp, int dissipator_id)
{
	crsMatrix * mat;

	mat = new crsMatrix(cd.Ns, cd.Ns);
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;
		Value[state_id_1].re = double(bit_at(cd.id_to_x[state_id_1], dissipator_id)) - double(bit_at(cd.id_to_x[state_id_1], dissipator_id + 1));
	}
	RowIndex[cd.Ns] = cd.Ns;

	return mat;
}

crsMatrix * create_A2_diss1_matrix(Model * m, ConfigData &cd, ConfigParam &cp, int dissipator_id)
{
	crsMatrix * mat;

	mat = new crsMatrix(cd.Ns, cd.Ns);
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;


	int tmp = 0;
	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		RowIndex[state_id_1] = state_id_1;

		tmp = bit_at(cd.id_to_x[state_id_1], dissipator_id) - bit_at(cd.id_to_x[state_id_1], dissipator_id + 1);

		if (tmp == 0)
		{
			Col[state_id_1] = state_id_1;
			Value[state_id_1].re = 0.0;
			Value[state_id_1].im = 0.0;
		}
		else
		{
			for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
			{
				if (cd.adjacement[cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2]])
				{
					vector<int> adjacency_bits = convert_int_to_vector_of_bits(cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2], cd.Nc);
					vector<int> hop;
					for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
					{
						if (adjacency_bits[cell_id])
						{
							hop.push_back(cell_id);
						}
					}

					for (int ad_cell_id = 0; ad_cell_id < hop.size(); ad_cell_id++)
					{
						hop[ad_cell_id] = (cd.Nc - 1) - hop[ad_cell_id];
					}

					if (hop[1] == dissipator_id)
					{
						if (bit_at(cd.id_to_x[state_id_1], dissipator_id))
						{
							Col[state_id_1] = state_id_2;
							Value[state_id_1].re = sin(-cp.dp);
							Value[state_id_1].im = -cos(-cp.dp);
						}
						else
						{
							Col[state_id_1] = state_id_2;
							Value[state_id_1].re = -sin(cp.dp);
							Value[state_id_1].im = cos(cp.dp);
						}
					}

					adjacency_bits.clear();
					hop.clear();
				}
			}
		}
	}
	RowIndex[cd.Ns] = cd.Ns;

	return mat;
}

void init_a1_a2_diss1(Model * m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;
	FMatrixs *Fs = m->Fs;

	int N_mat = (N + 1) * (N + 1) - 1;

	crsMatrix * result_a_matrix = NULL;

	crsMatrix * A1 = create_A1_diss1_matrix(m, cd, cp, 0);
	crsMatrix * A2 = create_A2_diss1_matrix(m, cd, cp, 0);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	int k = 0;

	crsMatrix * res;
	int cnt;
	a1_mat->RowIndex[0] = 0;
	for(int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A1, *(Fs->F[i+1]), *res);
		cnt = trace_struct(*res);
		if(cnt > 0)
		{
			a1_mat->Value[k] = trace(*res);
			a1_mat->Col[k] = i;
			k++;
		}
		a1_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a1_mat->NZ = k;

	k = 0;
	a2_mat->RowIndex[0] = 0;
	for(int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A2, *(Fs->F[i+1]), *res);
		cnt = trace_struct(*res);
		if(cnt > 0)
		{
			a2_mat->Value[k] = trace(*res);
			a2_mat->Col[k] = i;
			k++;
		}
		a2_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a2_mat->NZ = k;

	vector<pair<int , dcomplex> > * a_std;
	a_std = new vector<pair<int , dcomplex> >[N_mat];

	for(int i = 0; i < N_mat; i++)
	{
		for(int j = 0; j < N_mat; j++)
		{
			dcomplex v1, v2, v3, v4, res1,res2, res;
			int c1, c2, c3, c4;
			c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
			c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
			c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
			c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];

			if((c1 + c3) * (c2 + c4) > 0)
			{

				v1.re = a1_mat->Value[a1_mat->RowIndex[i]].re * c1;
				v3.re = a2_mat->Value[a2_mat->RowIndex[i]].re * c3;
				v2.re = a1_mat->Value[a1_mat->RowIndex[j]].re * c2;
				v4.re = a2_mat->Value[a2_mat->RowIndex[j]].re * c4;

				v1.im = a1_mat->Value[a1_mat->RowIndex[i]].im * c1;
				v3.im = a2_mat->Value[a2_mat->RowIndex[i]].im * c3;
				v2.im = a1_mat->Value[a1_mat->RowIndex[j]].im * c2;
				v4.im = a2_mat->Value[a2_mat->RowIndex[j]].im * c4;

				res1.re = v1.re + v3.im;
				res1.im = v1.im - v3.re;

				res2.re = v2.re + v4.im;
				res2.im = v2.im - v4.re;

				res2.im = -res2.im;

				res.re = res1.re * res2.re - res1.im * res2.im;
				res.im = res1.re * res2.im + res1.im * res2.re;
				a_std[i].push_back(make_pair(j, res));
			}
		}
	}

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a_mat_tmp_1 = NULL;
	a_mat_tmp_1 = stdToCrs(a_std, N_mat);
	result_a_matrix = new crsMatrix(*a_mat_tmp_1);

	delete a_mat_tmp_1;

	int diss_num = cd.Nc - 1;
	if (cp.bc == 1)
	{
		diss_num = cd.Nc - 1;
	}

	if (N > 1)
	{
		for (int dissId = 1; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss1_matrix(m, cd, cp, dissId);
			crsMatrix * A2 = create_A2_diss1_matrix(m, cd, cp, dissId);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			int k = 0;

			crsMatrix * res;
			int cnt;
			a1_mat->RowIndex[0] = 0;
			for(int i = 0; i < N_mat; i++)
			{
				res = new crsMatrix;
				SparseMKLMult(*A1, *(Fs->F[i+1]), *res);
				cnt = trace_struct(*res);
				if(cnt > 0)
				{
					a1_mat->Value[k] = trace(*res);
					a1_mat->Col[k] = i;
					k++;
				}
				a1_mat->RowIndex[i + 1] = k;
				delete res;
			}
			a1_mat->NZ = k;

			k = 0;
			a2_mat->RowIndex[0] = 0;
			for(int i = 0; i < N_mat; i++)
			{
				res = new crsMatrix;
				SparseMKLMult(*A2, *(Fs->F[i+1]), *res);
				cnt = trace_struct(*res);
				if(cnt > 0)
				{
					a2_mat->Value[k] = trace(*res);
					a2_mat->Col[k] = i;
					k++;
				}
				a2_mat->RowIndex[i + 1] = k;
				delete res;
			}
			a2_mat->NZ = k;

			vector<pair<int , dcomplex> > * a_std;
			a_std = new vector<pair<int , dcomplex> >[N_mat];

			for(int i = 0; i < N_mat; i++)
			{
				for(int j = 0; j < N_mat; j++)
				{
					dcomplex v1, v2, v3, v4, res1,res2, res;
					int c1, c2, c3, c4;
					c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
					c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
					c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
					c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];

					if((c1 + c3) * (c2 + c4) > 0)
					{

						v1.re = a1_mat->Value[a1_mat->RowIndex[i]].re * c1;
						v3.re = a2_mat->Value[a2_mat->RowIndex[i]].re * c3;
						v2.re = a1_mat->Value[a1_mat->RowIndex[j]].re * c2;
						v4.re = a2_mat->Value[a2_mat->RowIndex[j]].re * c4;

						v1.im = a1_mat->Value[a1_mat->RowIndex[i]].im * c1;
						v3.im = a2_mat->Value[a2_mat->RowIndex[i]].im * c3;
						v2.im = a1_mat->Value[a1_mat->RowIndex[j]].im * c2;
						v4.im = a2_mat->Value[a2_mat->RowIndex[j]].im * c4;

						res1.re = v1.re + v3.im;
						res1.im = v1.im - v3.re;

						res2.re = v2.re + v4.im;
						res2.im = v2.im - v4.re;

						res2.im = -res2.im;

						res.re = res1.re * res2.re - res1.im * res2.im;
						res.im = res1.re * res2.im + res1.im * res2.re;
						a_std[i].push_back(make_pair(j, res));
					}
				}
			}

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a_mat_tmp_1 = NULL;
			a_mat_tmp_1 = stdToCrs(a_std, N_mat);

			crsMatrix * a_mat_tmp_2 = new crsMatrix(*result_a_matrix);

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp_1, beta, *a_mat_tmp_2, *result_a_matrix);

			delete a_mat_tmp_1;
			delete a_mat_tmp_2;
		}
	}

	m->a_mat = new crsMatrix(*result_a_matrix);

	delete result_a_matrix;
}