#include "calcODE.h"
#include "genMatrix.h"
#include <math.h>
#include <stdlib.h>

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res)
{
	int i, j, s, f;
	for(i = 0; i < mat->N; i++)
	{
		s = mat->RowIndex[i];
		f = mat->RowIndex[i + 1];
		res[i].re = 0.0;
		res[i].im = 0.0;
		for(j = s; j < f; j++)
		{
			dcomplex v1 = mat->Value[j];
			dcomplex v2 = x[mat->Col[j]];
			res[i].re += v1.re * v2.re;// - v1.im * v2.im;
			res[i].im =0.0;//+= v1.re * v2.im + v1.im * v2.re;
		}
	}
}

void calcVectValue(double t, double h, 
				   Model * m, dcomplex *x, dcomplex * res, 
				   dcomplex * tmp1, dcomplex * tmp2)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * Gs = m->Gs;
	crsMatrix * QEs = m->QEs;
	dcomplex  * Ks = m->Ks;
	//double T = m->conf.T;
	//double A0 = m->conf.A0;
	//double w = m->conf.w;
	double T = 0.0;
	double A0 = 0.0;
	double w = 0.0;

	multMatVec(Gs, x, tmp1);
	multMatVec(QEs, x, tmp2);

	for(i = 0; i < N_mat; i++)
	{
		//res[i].re = (tmp1[i].re + A0 * sin(w * t) * tmp2[i].re - Ks[i].re) * h;
		/*if(sin(w * t) < 0.0)
		{
			res[i].re = (tmp1[i].re + A0 * (-1.0) * tmp2[i].re - Ks[i].re) * h;
		}
		else
		{
			res[i].re = (tmp1[i].re + A0 * (+1.0) * tmp2[i].re - Ks[i].re) * h;
		}*/
		//res[i].im = (tmp1[i].im + A0 * cos(w * t) * tmp2[i].im - Ks[i].im) * h;
		res[i].re = (tmp1[i].re + A0 * cos(w * t) * tmp2[i].re - Ks[i].re) * h;
		res[i].im = 0.0;//(tmp1[i].im + A0 * sin(w * t + PI) * tmp2[i].im - Ks[i].im) * h;
	}
}

dcomplex calcDiffIter(Model *m)
{
	dcomplex max_diff, diff;
	max_diff.re = 0.0;
	max_diff.im = 0.0;
	for(int i = 0; i < m->N_mat; i++)
	{
		diff.re = fabs(m->prevRhoF[i].re - m->RhoF[i].re);
		diff.im = fabs(m->prevRhoF[i].im - m->RhoF[i].im);
		if(max_diff.re < diff.re)max_diff.re = diff.re;
		if(max_diff.im < diff.im)max_diff.im = diff.im;
	}

	return max_diff;
}

void calcODE(Model *m, double h, int cntItr, double t)
{
	double time;
	int itr, i;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;

	for(i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1   = new dcomplex[N_mat];
	dcomplex * k2   = new dcomplex[N_mat];
	dcomplex * k3   = new dcomplex[N_mat];
	dcomplex * k4   = new dcomplex[N_mat];
	dcomplex * val  = new dcomplex[N_mat];
	dcomplex * tmp1 = new dcomplex[N_mat];
	dcomplex * tmp2 = new dcomplex[N_mat];
	dcomplex * tmp3 = new dcomplex[N_mat];

	for(itr = 0; itr < cntItr; itr ++)
	{
		time = t + h * itr;
		calcVectValue(time, h, m, RhoF, k1, tmp1, tmp2);
		for(i = 0; i < N_mat; i++)
		{
			val[i].re = RhoF[i].re + k1[i].re / 2.0;
			val[i].im = RhoF[i].im + k1[i].im / 2.0;
		}
		calcVectValue(time + h / 2.0, h, m, val, k2, tmp1, tmp2);
		for(i = 0; i < N_mat; i++)
		{
			val[i].re = RhoF[i].re + k2[i].re / 2.0;
			val[i].im = RhoF[i].im + k2[i].im / 2.0;
		}
		calcVectValue(time + h / 2.0, h, m, val, k3, tmp1, tmp2);
		for(i = 0; i < N_mat; i++)
		{
			val[i].re = RhoF[i].re + k3[i].re;
			val[i].im = RhoF[i].im + k3[i].im;
		}
		calcVectValue(time + h, h, m, val, k4, tmp1, tmp2);

		for(i = 0; i < N_mat; i++)
		{
			RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
			RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
		}
	}

	delete [] k1  ;
	delete [] k2  ;
	delete [] k3  ;
	delete [] k4  ;
	delete [] val ;
	delete [] tmp1;
	delete [] tmp2;
	delete [] tmp3;
}

dcomplex multMatCRS_tr(dcomplex *a, crsMatrix *b)
{
	int N = b->N;

	dcomplex c;
	c.re = 0.0;
	c.im = 0.0;

	crsMatrix *tB = new crsMatrix(*b);
	Transpose(*b, *tB,false);

	int i, j, jj, k, s, f;
	for(i = 0; i < N; i++)
	{
		s = tB->RowIndex[i];
		f = tB->RowIndex[i+1];
		for(k = s; k < f; k++)
		{
			j = tB->Col[k];
			c.re += a[i * N + j].re * tB->Value[k].re - a[i * N + j].im * tB->Value[k].im;
			c.im += a[i * N + j].re * tB->Value[k].im + a[i * N + j].im * tB->Value[k].re;
		}
	}

	delete tB;

	return c;
}

void initRhoODE(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	FMatrixs  * Fs = m->Fs;
	dcomplex  * RhoF = m->RhoF;
	dcomplex * psi = new dcomplex[(N + 1)*(N + 1)];
	createInitialMatrix(N + 1, psi);
	//createHermitianMatrix(N + 1, psi);

	for(int i = 0; i < N_mat; i++)
	{
		RhoF[i] = multMatCRS_tr(psi, Fs->F[i+1]);
	}

	delete[] psi;
}
