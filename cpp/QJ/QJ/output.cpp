#include "header.h"
#include "output.h"





MKL_Complex16 result_ksi_one_track(MKL_Complex16 * c, MKL_Complex16 * a, int N)
{
	MKL_Complex16 ZERO={0,0}, ONE={1,0}, I={0,1};
	MKL_Complex16 * tmp = new MKL_Complex16[N];
	MKL_Complex16 * c_tmp = new MKL_Complex16[N];

	memset(c_tmp, 0, N * sizeof(MKL_Complex16));
	double norm = 0;

	norm=sqrt(norm_vector2(c, N));
	for(int j = 0; j < N; j++)
	{
		c_tmp[j].real += c[j].real/norm;
		c_tmp[j].imag += c[j].imag/norm;
	}

	MKL_Complex16 x = {0,0};

	cblas_zgemv(CblasRowMajor, CblasNoTrans, N, N, &ONE, a, N, c_tmp, 1,  &ZERO, tmp, 1);

	for(int i = 0; i < N; i++)
	{
		x.real += c_tmp[i].real * tmp[i].real + c_tmp[i].imag * tmp[i].imag;
		x.imag += c_tmp[i].real * tmp[i].imag - c_tmp[i].imag * tmp[i].real;
	}

	delete(tmp);
	delete(c_tmp);

	return x;
}