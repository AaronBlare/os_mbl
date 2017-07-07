#pragma once
#include "config.h"
#include "data.h"
#include "mkl_vsl.h"

int SparseMKLMult(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLMultOne(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);

int SparseMKLAdd(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddOne(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);

void toOneBase(crsMatrix &A);
void toZeroBase(crsMatrix &A);
void Transpose(crsMatrix &Mat, crsMatrix &TMat, bool conj = true);
dcomplex trace(crsMatrix &A);
int trace_struct(crsMatrix &A);
void printMatrix(crsMatrix *A);
void printMatrixVal(crsMatrix *A);
void printVectorVal(dcomplex *A, int N);
void saveAngleMatrixVal(char* file, crsMatrix *A);
void saveAbsMatrixVal(char* file, crsMatrix *A);
void AbsMatrixDiagVal(crsMatrix *A, double * diag);
void saveVectorVal(char* file, dcomplex *vec, int N, int M);
void saveMatrix(char* file, crsMatrix *A);

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

void initFs(FMatrixs *Fs, int N);
void outFs(FMatrixs *Fs);
crsMatrix * createFeyeType(int N);
crsMatrix * createFPairTypeRe(int N, int i, int j);
crsMatrix * createFPairTypeIm(int N, int i, int j);
crsMatrix * createLastType(int N, int i);

void init_hE_vector(Model * m, ConfigData &cd, ConfigParam &cp);
void init_hU_vector(Model * m, ConfigData &cd, ConfigParam &cp);
void init_hJ_vector(Model * m, ConfigData &cd, ConfigParam &cp);
void init_H0(Model * m);

crsMatrix * stdToCrs(vector<pair<int, dcomplex> > * mat, int N);

void init_a1_a2_diss1(Model * m, ConfigData &cd, ConfigParam &cp);

struct Tensor_Coordinates
{
	MKL_Complex16 * data;
	unsigned int * coord1;
	unsigned int * coord2;
	unsigned int * coord3;
	unsigned int * hash;
	unsigned int k;
};

struct Tensor_Coordinates_1
{
	MKL_Complex16 * data;
	unsigned int coord1;
	unsigned int * coord2;
	unsigned int * coord3;
	unsigned int N;
};

void sort_matrix(Tensor_Coordinates * matrix);
void fijk_coord(Tensor_Coordinates * f_ijk, int N);
void dijk_coord(Tensor_Coordinates * d_ijk, int N);
void free_matrix(Tensor_Coordinates * matrix);
void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat);

void init_f_d(Model *m);
void init_f_d_valentin(Model *m);

void transpFs(Model *m);

void calcQEs(Model * m);
void calcQUs(Model * m);
void calcQJs(Model * m);

void calcKs(Model *m, ConfigData &cd, ConfigParam &cp);

void calcRs(Model *m, ConfigData &cd, ConfigParam &cp);
crsMatrix* calcSubRs(Model *m, int start, int finish);

void calcGs(Model *m);

int createHermitianMatrix(int n, dcomplex *a);
int createInitialMatrix(int n, dcomplex *a);
int genNormalDistributedElemets(int n1, int n2, double * re, double * im);

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res);
void initRhoODE(Model *m, ConfigData &cd, ConfigParam &cp);
void calcODE(Model *m, IntData &int_data, ConfigData &cd, ConfigParam &cp);
dcomplex calcDiffIter(Model *m);

void linSolv(Model *m);
void linSolvCheck(Model *m);
void linSolvReal(Model *m);

void calcRho(Model *m);
void clearRho(Model *m);

void characteristics(Model *m, ConfigData &cd, ConfigParam &cp, bool append);
