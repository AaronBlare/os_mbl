#include "CalcRs.h"
#include "data.h"
#include "config.h"

#include <omp.h>
crsMatrix* calcSubRs(Model *m, int start, int finish);

void calcRs(Model *m, ConfigData &cd, ConfigParam &cp)
{
	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	int N_mat = m->N_mat;

	crsMatrix * Rs_tmp;
	crsMatrix ** RsTh;
	crsMatrix * Rs = new crsMatrix(m->N_mat, 1);

	Rs->Col[0] = 1;
	Rs->RowIndex[0] = 1;
	for(int i = 1; i <= N_mat; i++)
	{
		Rs->RowIndex[i] = 1;
	}

	int nThread = 1;
	//  omp_set_num_threads(3);
#pragma omp parallel
	{
#pragma omp single
		nThread = omp_get_num_threads();
	}
	printf("nThread %d \n", nThread);

	int SubTask = 1;
	int CntTask = nThread * SubTask;

	RsTh = new crsMatrix *[CntTask];

	crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_H_mat = m->f_H_mat;
	crsMatrix ** d_mat = m->d_mat;

	for(int i = 0; i < N_mat; i++)
	{
		toOneBase(*(f_mat[i]));
		toOneBase(*(d_mat[i]));
		toOneBase(*(f_H_mat[i]));
	}

	int size = N_mat / CntTask;
	//  int *start = new int[CntTask];
	//  int *finish = new int[CntTask];

	//  start[0] = 0;
	//  finish[0] = size;
	//  if((N_mat % CntTask) != 0) finish[0]++;
	printf("N_mat %d \n", N_mat);
	//  printf("%d %d \n", start[0], finish[0]);
	//  for(int i = 1; i < CntTask; i++)
	//  {
	//    start[i] = finish[i - 1];
	//    finish[i] = start[i] + size;
	//    if((N_mat % CntTask) >  i)finish[i]++;
	//    printf("%d %d \n", start[i], finish[i]);
	//  }


	int id = 0;
	int portion = 100;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int localTskID;
		while (true)
		{
#pragma omp critical
			{
				localTskID = id;
				id++;
			}
			int start = localTskID * portion;
			int finish = (localTskID + 1) * portion;
			if (start > N_mat) break;
			if (finish > N_mat) {
				finish = N_mat;
				printf("%d %d \n", start, finish);
			}
			RsTh[tid] = calcSubRs(m, start, finish);
#pragma omp critical
			{
				Rs_tmp = new crsMatrix;
				SparseMKLAddOne(*Rs, sum, *RsTh[tid], *Rs_tmp);
				delete RsTh[tid];
				delete Rs;
				Rs = Rs_tmp;
			}
		}

	}
	//delete [] start;
	//delete [] finish;

	//  for(int i = 0; i < CntTask; i++)
	//  {
	//      Rs_tmp = new crsMatrix;
	//      SparseMKLAddOne(*Rs, sum, *RsTh[i], *Rs_tmp);
	//      delete RsTh[i];
	//      delete Rs;
	//      Rs = Rs_tmp;
	//  }
	delete[] RsTh;

	//Rs = calcSubRs(m, 0, N_mat);
	toZeroBase(*Rs);
	for(int i = 0; i < N_mat; i++)
	{
		toZeroBase(*(f_mat[i]));
		toZeroBase(*(d_mat[i]));
		toZeroBase(*(f_H_mat[i]));
	}
	//  delete FkT;
	//  delete FiT;

	double mm = -cp.g / 4.0;
	for(int i = 0; i < Rs->NZ; i++)
	{
		Rs->Value[i].re *= mm;
		Rs->Value[i].im *= mm;
	}
	m->Rs = Rs;
}


crsMatrix* calcSubRs(Model *m, int start, int finish)
{
	int N_mat = m->N_mat;
	dcomplex sum_i;
	sum_i.re = 0.0;
	sum_i.im = 1.0;

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	crsMatrix * tmp;
	crsMatrix * MatSDi, * MatSDk;
	crsMatrix * M1, *M2, *Rsum;

	crsMatrix * Rs_tmp;
	crsMatrix * Rs = new crsMatrix(m->N_mat, 1);

	Rs->Col[0] = 1;
	Rs->RowIndex[0] = 1;
	for(int i = 1; i <= N_mat; i++)
	{
		Rs->RowIndex[i] = 1;
	}

	crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_H_mat = m->f_H_mat;
	crsMatrix ** d_mat = m->d_mat;
	crsMatrix *  As    = m->a_mat;

	M1 = new crsMatrix;
	M2 = new crsMatrix;
	MatSDi = new crsMatrix;
	MatSDk = new crsMatrix;
	Rsum = new crsMatrix;

	//  for(int i = 0; i < N_mat; i++)
	//  {
	//    toOneBase(*(f_mat[i]));
	//    toOneBase(*(d_mat[i]));
	//    toOneBase(*(f_H_mat[i]));
	//  }

	for(int i = start; i < finish; i++)
	{
		SparseMKLAddOne(*(f_mat[i]), sum_i, *(d_mat[i]), *MatSDi, true);

		int st = As->RowIndex[i];
		int fn = As->RowIndex[i + 1];

		for(int j = st; j < fn; j++)
		{
			int k = As->Col[j];
			SparseMKLAddOne(*(f_mat[k]), sum_i, *(d_mat[k]), *MatSDk, true);

			for(int conj_i = 0; conj_i < MatSDk->NZ; conj_i++)
				MatSDk->Value[conj_i].im = -MatSDk->Value[conj_i].im;

			SparseMKLMultOne(*(f_H_mat[k]), *MatSDi, *M1, true);
			SparseMKLMultOne(*(f_H_mat[i]), *MatSDk, *M2, true);

			SparseMKLAddOne(*M1, sum, *M2, *Rsum, true);

			dcomplex val1 = As->Value[j];
			for(int ii = 0; ii < Rsum->NZ; ii ++)
			{
				dcomplex val2 = Rsum->Value[ii];
				Rsum->Value[ii].re = val1.re *  val2.re - val1.im *  val2.im;
				Rsum->Value[ii].im = val1.re *  val2.im + val1.im *  val2.re;
			}

			Rs_tmp = new crsMatrix;
			SparseMKLAddOne(*Rs, sum, *Rsum, *Rs_tmp);
			delete Rs;
			Rs = Rs_tmp;
		}
	}

	//  toZeroBase(*Rs);
	//  for(int i = 0; i < N_mat; i++)
	//  {
	//    toZeroBase(*(f_mat[i]));
	//    toZeroBase(*(d_mat[i]));
	//    toZeroBase(*(f_H_mat[i]));
	//  }
	//  delete FkT;
	//  delete FiT;

	delete Rsum;

	delete MatSDi;
	delete MatSDk;

	delete M1;
	delete M2;

	//double mm = -m->conf.g / 4.0 / m->N;
	//for(int i = 0; i < Rs->NZ; i++)
	//{
	//  Rs->Value[i].re *= mm;
	//  Rs->Value[i].im *= mm;
	//}
	return Rs;
}