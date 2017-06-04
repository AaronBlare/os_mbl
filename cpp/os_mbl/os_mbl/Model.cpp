#include "Model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Model * createModel(int N, ConfigParam conf)
{
	int i;
	Model * model = new Model;

	model->N = N;
	model->N_mat = (N + 1) * (N + 1) - 1;
	model->conf = conf;

	model->Fs = new FMatrixs;
	createFMatrixs(model->Fs, N);

	model->hE = new dcomplex[model->N_mat];
	model->hU = new dcomplex[model->N_mat];
	model->hJ = new dcomplex[model->N_mat];

	model->HE = NULL;
	model->HU = NULL;
	model->HJ = NULL;

	model->H0 = NULL;

	model->f_mat = NULL;
	model->f_H_mat = NULL;
	model->d_mat = NULL;

	model->a_mat = NULL;

	model->QEs = NULL;
	model->QUs = NULL;
	model->QJs = NULL;

	model->Ks = new dcomplex[model->N_mat];
	model->Rs = NULL;
	model->Gs = NULL;

	model->prevRhoF = new dcomplex[model->N_mat];
	model->RhoF = new dcomplex[model->N_mat];
	memset(model->RhoF, 0, sizeof(dcomplex) * model->N_mat);
	memset(model->prevRhoF, 0, sizeof(dcomplex) * model->N_mat);
	for (i = 0; i < model->N_mat; i++)
	{
		model->RhoF[i].re = (double)rand() / (double)RAND_MAX;
		model->RhoF[i].im = 0.0;
	}

	model->Rho = NULL;

	return model;
}

void createFMatrixs(FMatrixs * Fs, int N)
{
	Fs->countF = (2 + N) * N + 1;
	Fs->F = new crsMatrix *[Fs->countF];
	for (int i = 0; i < Fs->countF; i++)
	{
		Fs->F[i] = NULL;
	}
}

//освобождение памяти из под модели
void freeModel(Model * model)
{
	freeFMatrixs(model->Fs);
	delete model->Fs;
	delete[] model->hE;
	delete[] model->hU;
	delete[] model->hJ;

	if (model->f_mat != NULL)
	{
		int N = model->N;
		for (int i = 0; i < (N + 1) * (N + 1) - 1; i++)
		{
			delete model->f_mat[i];
		}
		delete[] model->f_mat;
	}

	if (model->f_H_mat != NULL)
	{
		int N = model->N;
		for (int i = 0; i < (N + 1) * (N + 1) - 1; i++)
		{
			delete model->f_H_mat[i];
		}
		delete[] model->f_H_mat;
	}

	if (model->d_mat != NULL)
	{
		int N = model->N;
		for (int i = 0; i < (N + 1) * (N + 1) - 1; i++)
		{
			delete model->d_mat[i];
		}
		delete[] model->d_mat;
	}

	if (model->a_mat != NULL)
	{
		delete model->a_mat;
	}

	if (model->QEs != NULL)
	{
		delete model->QEs;
	}

	if (model->QUs != NULL)
	{
		delete model->QUs;
	}

	if (model->QJs != NULL)
	{
		delete model->QJs;
	}

	if (model->HE != NULL)
	{
		delete model->HE;
	}

	if (model->HU != NULL)
	{
		delete model->HU;
	}

	if (model->HJ != NULL)
	{
		delete model->HJ;
	}

	if (model->H0 != NULL)
	{
		delete model->H0;
	}

	if (model->Ks != NULL)
	{
		delete[]model->Ks;
	}
	if (model->Rs != NULL)
	{
		delete model->Rs;
	}

	if (model->Gs != NULL)
	{
		delete model->Gs;
	}

	if (model->RhoF != NULL)
	{
		delete[] model->RhoF;
	}

	if (model->prevRhoF != NULL)
	{
		delete[] model->prevRhoF;
	}

	if (model->Rho != NULL)
	{
		delete model->Rho;
	}
}

void freeFMatrixs(FMatrixs * Fs)
{
	for (int i = 0; i < Fs->countF; i++)
	{
		if (Fs->F[i] != NULL)
		{
			delete Fs->F[i];
		}
	}
	delete[] Fs->F;
}