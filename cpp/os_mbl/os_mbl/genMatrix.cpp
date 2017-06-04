#include <stdlib.h>
#include "mkl_vsl.h"
#include "Model.h"
#include <math.h>

#define SEED    777
#define BRNG    VSL_BRNG_MCG31
#define METHOD  0

int allocComplexMatrix(int n, double ** re, double ** im)
{
	double * mem;
	mem = (double *) malloc(sizeof(double) * n * n);
	if(mem == NULL)
	{
		return -1;
	}
	(*re) = mem;

	mem = (double *) malloc(sizeof(double) * n * n);
	if(mem == NULL)
	{
		return -1;
	}
	(*im) = mem;

	return 0;
}

int freeComplexMatrix(double * re, double * im)
{
	free(re);
	free(im);
	return 0;
}
int genNormalDistributedElemets(int n1, int n2, double * re, double * im);

int genNormalDistributedElemets(int n, double * re, double * im)
{
	return genNormalDistributedElemets(n, n, re, im);
}
int genNormalDistributedElemets(int n1, int n2, double * re, double * im)
{
	VSLStreamStatePtr stream;
	int i, errcode;
	double a=0.0,sigma=1.0;

	/***** Initialize *****/
	errcode = vslNewStream( &stream, BRNG,  SEED );
	if(errcode != VSL_STATUS_OK) return errcode;

	errcode = vdRngGaussian( METHOD, stream, n1 * n2, re, a, sigma );
	if(errcode != VSL_STATUS_OK) return errcode;

	errcode = vdRngGaussian( METHOD, stream, n1 * n2, im, a, sigma );
	if(errcode != VSL_STATUS_OK) return errcode;

	/***** Deinitialize *****/
	errcode = vslDeleteStream( &stream );
	if(errcode != VSL_STATUS_OK) return errcode;

	return 0;
}

int createHermitianMatrix(int n, dcomplex *a)
{
	double *re1, *im1, *re2, *im2;
	int errcode, i, j;
	errcode = allocComplexMatrix(n, &re1, &im1);
	if(errcode != 0) return -3;

	errcode = allocComplexMatrix(n, &re2, &im2);
	if(errcode != 0) return -4;

	genNormalDistributedElemets(n, re1, im1);
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			re2[i*n + j] =  re1[j*n + i];
			im2[i*n + j] = -im1[j*n + i];
		}
	}

	double sre = 0.0, sim=0.0;

	double two = 2.0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			re1[i*n + j] = (re1[i*n + j] + re2[i*n + j]) / two;
			im1[i*n + j] = (im1[i*n + j] + im2[i*n + j]) / two;
			sre += fabs(re1[i*n + j]);
			sim += fabs(im1[i*n + j]);
		}
	}

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			re1[i*n + j] /= sre;
			im1[i*n + j] /= sim;
		}
	}

	freeComplexMatrix(re2, im2);  

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			a[i*n + j].re = re1[i*n + j];
			a[i*n + j].im = im1[i*n + j];
		}
	}

	freeComplexMatrix(re1, im1);

	return 0;
}


int createInitialMatrix(int n, dcomplex *a)
{
	double *re, *im;
	int errcode, i, j;
	//double sre = 0.0, sim=0.0;
	double sum = 0.0;

	errcode = allocComplexMatrix(n, &re, &im);
	if(errcode != 0) return -3;

	genNormalDistributedElemets(n, 1, re, im);
	for(i = 0; i < n; i++)
	{
		//sre += abs(re[i]);
		//sim += abs(im[i]);
		sum += re[i] * re[i] + im[i] * im[i];
	}

	for(i = 0; i < n; i++)
	{
		//re[i] /= sre;
		//im[i] /= sim;
		//re[i] /= sqrt(sum);
		//im[i] /= sqrt(sum);
		re[i] = sqrt(1.0 / n / 2);
		im[i] = sqrt(1.0 / n / 2);
	}
	sum = 0.0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			a[i*n + j].re = re[i]*re[j] + im[i]*im[j];
			a[i*n + j].im = re[i]*im[j] - im[i]*re[j];
		}
		sum += a[i*n + i].re;
	}

	printf("Init Tr(Rho): %lf\n", sum);

	freeComplexMatrix(re, im);

	return 0;

}
