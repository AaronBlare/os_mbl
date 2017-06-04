#include "CalcQs.h"

void calcQEs(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * hE = m->hE;

	//crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * QEs = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	QEs->Col[0] = 0;
	QEs->RowIndex[0] = 0;
	for(int i = 1; i <= N_mat; i++)
	{
		QEs->RowIndex[i] = 0;
	}

	for(int i = 0; i < N_mat; i++)
	{
		if((hE[i].re != 0.0) || (hE[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*QEs, hE[i], *(f_mat[i]), *resSum);
			delete QEs;
			QEs = resSum;
		}
	}

	m->QEs = QEs;
}

void calcQUs(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * hU = m->hU;

	//crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * QUs = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	QUs->Col[0] = 0;
	QUs->RowIndex[0] = 0;
	for(int i = 1; i <= N_mat; i++)
	{
		QUs->RowIndex[i] = 0;
	}

	for(int i = 0; i < N_mat; i++)
	{
		if((hU[i].re != 0.0) || (hU[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*QUs, hU[i], *(f_mat[i]), *resSum);
			delete QUs;
			QUs = resSum;
		}
	}

	m->QUs = QUs;
}

void calcQJs(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * hJ = m->hJ;

	//crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * QJs = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	QJs->Col[0] = 0;
	QJs->RowIndex[0] = 0;
	for(int i = 1; i <= N_mat; i++)
	{
		QJs->RowIndex[i] = 0;
	}

	for(int i = 0; i < N_mat; i++)
	{
		if((hJ[i].re != 0.0) || (hJ[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*QJs, hJ[i], *(f_mat[i]), *resSum);
			delete QJs;
			QJs = resSum;
		}
	}

	m->QJs = QJs;
}
