#ifndef __MODEL__
#define __MODEL__

#pragma once
#include "Matrix_op.h"
#include "config.h"
#include "data.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

struct FMatrixs
{
	crsMatrix **F;
	int countF;
};

struct Model
{
	int N;
	int N_mat;
	ConfigParam conf;
	FMatrixs *Fs;

	dcomplex *hE;
	dcomplex *hU;
	dcomplex *hJ;

	crsMatrix * HE;
	crsMatrix * HU;
	crsMatrix * HJ;

	crsMatrix * H0;

	dcomplex *he;

	crsMatrix ** f_mat;
	crsMatrix ** f_H_mat;
	crsMatrix ** d_mat;

	crsMatrix * a_mat;

	crsMatrix * QEs;
	crsMatrix * QUs;
	crsMatrix * QJs;

	crsMatrix * Qs;
	dcomplex  * Ks;
	crsMatrix * Rs;
	crsMatrix * Gs;

	dcomplex  * prevRhoF;
	dcomplex  * RhoF;

	crsMatrix * Rho;
};

Model * createModel(int N, ConfigParam conf);
void freeModel(Model * model);

void createFMatrixs(FMatrixs * Fs, int N);
void freeFMatrixs(FMatrixs * Fs);

#endif