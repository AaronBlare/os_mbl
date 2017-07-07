#include "f_basis.h"
#include <string.h>


#define IND(i, j, k) (i)*(N * N - 1)*(N * N - 1) + (j)*(N * N - 1) + (k)
#define IndS(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2)
#define IndJ(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2 + 1)
#define IndD(l)    (N * (N-1) + l - 1)

#define SEED    777
#define BRNG    VSL_BRNG_MCG31
#define METHOD  0

int SparseMKLMult(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
		B.RowIndex[j]++;
	}

	// Используется функция, вычисляющая C = op(A) * B
	char trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

					  // Хитрый параметр, влияющий на то, как будет выделяться память
					  // request = 0: память для результирующей матрицы д.б. выделена заранее
					  // Если мы не знаем, сколько памяти необходимо для хранения результата,
					  // необходимо:
					  // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
					  // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
					  //                                                         последний элемент
					  // 3) выделить память для массивов c и jc 
					  //    (кол-во элементов = ic[Кол-во строк]-1)
					  // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;

	request = 1;
	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C

		C.Value = 0;
		C.Col = 0;
	}
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]--;
	for (i = 0; i < C.NZ; i++)
		C.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
		B.RowIndex[j]--;
		C.RowIndex[j]--;
	}

	return 0;
}
int SparseMKLMultOne(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;

	// Используется функция, вычисляющая C = op(A) * B
	char trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

					  // Хитрый параметр, влияющий на то, как будет выделяться память
					  // request = 0: память для результирующей матрицы д.б. выделена заранее
					  // Если мы не знаем, сколько памяти необходимо для хранения результата,
					  // необходимо:
					  // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
					  // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
					  //                                                         последний элемент
					  // 3) выделить память для массивов c и jc 
					  //    (кол-во элементов = ic[Кол-во строк]-1)
					  // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;

	request = 1;
	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C

		C.Value = 0;
		C.Col = 0;
	}
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С

	return 0;
}

int SparseMKLAdd(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
		B.RowIndex[j]++;
	}

	// Используется функция, вычисляющая C = op(A) * B
	char trans;

	trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

				 // Хитрый параметр, влияющий на то, как будет выделяться память
				 // request = 0: память для результирующей матрицы д.б. выделена заранее
				 // Если мы не знаем, сколько памяти необходимо для хранения результата,
				 // необходимо:
				 // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
				 // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
				 //                                                         последний элемент
				 // 3) выделить память для массивов c и jc 
				 //    (кол-во элементов = ic[Кол-во строк]-1)
				 // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;
	request = 1;

	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C
		C.Value = 0;
		C.Col = 0;
	}

	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]--;
	for (i = 0; i < C.NZ; i++)
		C.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
		B.RowIndex[j]--;
		C.RowIndex[j]--;
	}

	return 0;
}
int SparseMKLAddT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
		B.RowIndex[j]++;
	}

	// Используется функция, вычисляющая C = op(A) * B
	char trans;

	trans = 'T'; // говорит о том, op(A) = A - не нужно транспонировать A

	int request;

	int sort = 8;

	int nzmax = -1;

	// Служебная информация
	int info;
	request = 1;

	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C
		C.Value = 0;
		C.Col = 0;
	}

	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]--;
	for (i = 0; i < C.NZ; i++)
		C.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
		B.RowIndex[j]--;
		C.RowIndex[j]--;
	}

	return 0;
}
int SparseMKLAddOne(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;

	// Используется функция, вычисляющая C = op(A) * B
	char trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

					  // Хитрый параметр, влияющий на то, как будет выделяться память
					  // request = 0: память для результирующей матрицы д.б. выделена заранее
					  // Если мы не знаем, сколько памяти необходимо для хранения результата,
					  // необходимо:
					  // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
					  // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
					  //                                                         последний элемент
					  // 3) выделить память для массивов c и jc 
					  //    (кол-во элементов = ic[Кол-во строк]-1)
					  // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;
	request = 1;

	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C
		C.Value = 0;
		C.Col = 0;
	}

	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С

	return 0;
}

void toOneBase(crsMatrix &A)
{
	int i, j, n = A.N;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
	}

}
void toZeroBase(crsMatrix &A)
{
	int i, j, n = A.N;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
	}
}
dcomplex trace(crsMatrix &A)
{
	dcomplex res;
	res.re = 0.0;
	res.im = 0.0;

	for (int i = 0; i < A.N; i++)
	{
		for (int k = A.RowIndex[i]; k < A.RowIndex[i + 1]; k++)
		{
			int j = A.Col[k];
			if (i == j)
			{
				res.re += A.Value[k].re;
				res.im += A.Value[k].im;
			}
		}
	}

	return res;
}
int trace_struct(crsMatrix &A)
{
	int res = 0;

	for (int i = 0; i < A.N; i++)
	{
		for (int k = A.RowIndex[i]; k < A.RowIndex[i + 1]; k++)
		{
			int j = A.Col[k];
			if (i == j)
			{
				res++;
			}
		}
	}

	return res;
}
void printMatrix(crsMatrix *A)
{
	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				printf("0 ");
				j++;
			}
			printf("1 ");
			j++;
		}
		while (j < A->N)
		{
			printf("0 ");
			j++;
		}
		printf("\n");
	}
	printf("\n");
}
void printMatrixVal(crsMatrix *A)
{
	printf("####################################\n");
	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				printf(" 0.0 ");
				j++;
			}
			printf("%3.5lf ", A->Value[k].re);
			j++;
		}
		while (j < A->N)
		{
			printf(" 0.0 ");
			j++;
		}
		printf("\n");
	}
	printf("\n");

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				printf(" 0.0 ");
				j++;
			}
			printf("%3.5lf ", A->Value[k].im);
			j++;
		}
		while (j < A->N)
		{
			printf(" 0.0 ");
			j++;
		}
		printf("\n");
	}
	printf("####################################\n");
	printf("\n");
}
void saveAbsMatrixVal(char* file, crsMatrix *A)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				fprintf(f, "0.00000000 ");
				j++;
			}
			fprintf(f, "%.8lf ", sqrt(A->Value[k].re * A->Value[k].re + A->Value[k].im * A->Value[k].im));
			j++;
		}
		while (j < A->N)
		{
			fprintf(f, "0.00000000 ");
			j++;
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
void AbsMatrixDiagVal(crsMatrix *A, double * diag)
{

	for (int i = 0; i < A->N; i++)
	{
		diag[i] = 0.0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			if (A->Col[k] == i)
			{
				diag[i] = sqrt(A->Value[k].re * A->Value[k].re + A->Value[k].im * A->Value[k].im);
			}
		}
	}
}
void saveAngleMatrixVal(char* file, crsMatrix *A)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				fprintf(f, "0.00000000 ");
				j++;
			}
			fprintf(f, "%.8lf ", atan2(A->Value[k].im, A->Value[k].re) * 180.0 / 3.14159265);
			j++;
		}
		while (j < A->N)
		{
			fprintf(f, "0.00000000 ");
			j++;
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
void saveVectorVal(char* file, dcomplex *vec, int N, int M)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fprintf(f, "%.10lf ", vec[i * M + j].re);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fprintf(f, "%.10lf ", vec[i * M + j].im);
		}
		fprintf(f, "\n");
	}

	fclose(f);
}
void printVectorVal(dcomplex *A, int N)
{
	printf("####################################\n");
	for (int i = 0; i < N; i++)
	{
		printf("%3.4lf ", A[i].re);
	}
	printf("\n");
	for (int i = 0; i < N; i++)
	{
		printf("%3.1lf ", A[i].im);
	}
	printf("\n");
	printf("####################################\n");
	printf("\n");
}
void Transpose(crsMatrix &Mat, crsMatrix &TMat, bool conj)
{
	int i, j, nz;
	int S;

	int n = Mat.N;
	int* column = Mat.Col;
	int* row = Mat.RowIndex;
	dcomplex* val = Mat.Value;

	nz = row[n];

	int* tColumn = TMat.Col;
	int* tRow = TMat.RowIndex;
	dcomplex* tVal = TMat.Value;

	memset(tRow, 0, (n + 1) * sizeof(int));
	for (i = 0; i < nz; i++)
		tRow[column[i] + 1]++;

	S = 0;
	for (i = 1; i <= n; i++)
	{
		int tmp = tRow[i];
		tRow[i] = S;
		S = S + tmp;
	}

	for (i = 0; i < n; i++)
	{
		int j1 = row[i];
		int j2 = row[i + 1];
		int Col = i; // Столбец в AT - строка в А
		for (j = j1; j < j2; j++)
		{
			dcomplex V = val[j];  // Значение
			int RIndex = column[j];  // Строка в AT
			int IIndex = tRow[RIndex + 1];
			tVal[IIndex] = V;
			tColumn[IIndex] = Col;
			tRow[RIndex + 1]++;
		}
	}
	if (conj)
	{
		for (i = 0; i < nz; i++)
		{
			tVal[i].im = -tVal[i].im;
		}
	}
}
void saveMatrix(char* file, crsMatrix *A)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{
			fprintf(f, "%d %d %.16lf %.16lf \n", i + 1, A->Col[k] + 1, A->Value[k].re, A->Value[k].im);
		}
	}
	fclose(f);
}

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

void outFs(FMatrixs *Fs)
{
	for (int i = 0; i < Fs->countF; i++)
	{
		printMatrixVal(Fs->F[i]);
	}
}
void initFs(FMatrixs *Fs, int N)
{
	int i, j, k;
	Fs->F[0] = createFeyeType(N + 1);
	k = 1;
	for (i = 0; i < N + 1; i++)
	{
		for (j = i + 1; j < N + 1; j++)
		{
			Fs->F[k] = createFPairTypeRe(N + 1, i, j); k++;
			Fs->F[k] = createFPairTypeIm(N + 1, i, j); k++;
		}
	}

	for (i = 0; i < N; i++)
	{
		if (k < Fs->countF)
		{
			Fs->F[k] = createLastType(N + 1, i); k++;
		}
		else
		{
			throw("error count calc (no mem)");
		}
	}

	if (k != Fs->countF)
	{
		throw("error count calc (countF > k)");
	}

	//outFs(Fs);
}
crsMatrix * createFeyeType(int N)
{
	crsMatrix * mat;
	mat = new crsMatrix(N, N);
	for (int i = 0; i < N; i++)
	{
		mat->Col[i] = i;
		mat->RowIndex[i] = i;
		mat->Value[i].re = 1.0;
	}
	mat->RowIndex[N] = N;

	return mat;
}
crsMatrix * createFPairTypeRe(int N, int i, int j)
{
	double val = 1.0 / sqrt(2.0);
	crsMatrix * mat;
	mat = new crsMatrix(N, 2);
	mat->Value[0].re = val;
	mat->Value[1].re = val;
	mat->Col[0] = j;
	mat->Col[1] = i;

	for (int ii = 0; ii < i + 1; ii++)
	{
		mat->RowIndex[ii] = 0;
	}
	for (int ii = i + 1; ii < j + 1; ii++)
	{
		mat->RowIndex[ii] = 1;
	}
	for (int ii = j + 1; ii <= N; ii++)
	{
		mat->RowIndex[ii] = 2;
	}

	return mat;
}
crsMatrix * createFPairTypeIm(int N, int i, int j)
{
	double val = -1.0 / sqrt(2.0);
	crsMatrix * mat;
	mat = new crsMatrix(N, 2);
	mat->Value[0].im = val;
	mat->Value[1].im = -val;
	mat->Col[0] = j;
	mat->Col[1] = i;

	for (int ii = 0; ii < i + 1; ii++)
	{
		mat->RowIndex[ii] = 0;
	}
	for (int ii = i + 1; ii < j + 1; ii++)
	{
		mat->RowIndex[ii] = 1;
	}
	for (int ii = j + 1; ii <= N; ii++)
	{
		mat->RowIndex[ii] = 2;
	}

	return mat;
}
crsMatrix * createLastType(int N, int i)
{
	crsMatrix * mat;
	mat = new crsMatrix(N, i + 2);
	int ii;
	double val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
	mat->RowIndex[0] = 0;
	for (ii = 0; ii <= i; ii++)
	{
		mat->Col[ii] = ii;
		mat->RowIndex[ii + 1] = mat->RowIndex[ii] + 1;
		mat->Value[ii].re = val;
	}
	mat->Col[ii] = ii;
	mat->RowIndex[ii + 1] = mat->RowIndex[ii] + 1;
	mat->Value[ii].re = -(i + 1) * val;
	ii++;

	for (; ii < N; ii++)
	{
		mat->RowIndex[ii + 1] = mat->RowIndex[ii];
	}

	return mat;
}

crsMatrix * create_HE_matrix(Model * m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * HE = new crsMatrix(N + 1, N + 1);
	int * RowIndex = HE->RowIndex;
	int * Col = HE->Col;
	dcomplex * Value = HE->Value;

	double * diag_array = new double[N + 1];

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
	for (i = 0; i < N + 1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for (i = 0; i < N + 1; i++)
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
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*HE, *(Fs->F[i + 1]), *res);
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

	double * hu_diag = new double[N + 1];

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		hu_diag[state_id] = cp.U * bit_count(cd.id_to_x[state_id] & (cd.id_to_x[state_id] << 1));
	}

	int i, j;
	RowIndex[0] = 0;
	for (i = 0; i < N + 1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for (i = 0; i < N + 1; i++)
	{
		Value[i].re = hu_diag[i];
	}

	delete[] hu_diag;

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
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*HU, *(Fs->F[i + 1]), *res);
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
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*HJ, *(Fs->F[i + 1]), *res);
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

crsMatrix * stdToCrs(vector<pair<int, dcomplex> > * mat, int N)
{
	crsMatrix * res = NULL;

	int NZ = 0;
	int i, j, k, Nl;

	for (i = 0; i < N; i++)
	{
		NZ += mat[i].size();
	}
	res = new crsMatrix(N, NZ);
	k = 0;
	res->RowIndex[0] = 0;
	for (i = 0; i < N; i++)
	{
		Nl = mat[i].size();
		for (j = 0; j < Nl; j++)
		{
			res->Col[k] = mat[i][j].first;
			res->Value[k] = mat[i][j].second;
			k++;
		}
		res->RowIndex[i + 1] = k;
	}

	return res;
}

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

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss1_matrix(m, cd, cp, 0);
	crsMatrix * A2 = create_A2_diss1_matrix(m, cd, cp, 0);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	int k = 0;

	crsMatrix * res;
	int cnt;
	a1_mat->RowIndex[0] = 0;
	for (int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A1, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
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
	for (int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A2, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
		{
			a2_mat->Value[k] = trace(*res);
			a2_mat->Col[k] = i;
			k++;
		}
		a2_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a2_mat->NZ = k;

	vector<pair<int, dcomplex> > * a_std;
	a_std = new vector<pair<int, dcomplex> >[N_mat];

	for (int i = 0; i < N_mat; i++)
	{
		for (int j = 0; j < N_mat; j++)
		{
			dcomplex v1, v2, v3, v4, res1, res2, res;
			int c1, c2, c3, c4;
			c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
			c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
			c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
			c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];

			if ((c1 + c3) * (c2 + c4) > 0)
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

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;

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
			for (int i = 0; i < N_mat; i++)
			{
				res = new crsMatrix;
				SparseMKLMult(*A1, *(Fs->F[i + 1]), *res);
				cnt = trace_struct(*res);
				if (cnt > 0)
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
			for (int i = 0; i < N_mat; i++)
			{
				res = new crsMatrix;
				SparseMKLMult(*A2, *(Fs->F[i + 1]), *res);
				cnt = trace_struct(*res);
				if (cnt > 0)
				{
					a2_mat->Value[k] = trace(*res);
					a2_mat->Col[k] = i;
					k++;
				}
				a2_mat->RowIndex[i + 1] = k;
				delete res;
			}
			a2_mat->NZ = k;

			vector<pair<int, dcomplex> > * a_std;
			a_std = new vector<pair<int, dcomplex> >[N_mat];

			for (int i = 0; i < N_mat; i++)
			{
				for (int j = 0; j < N_mat; j++)
				{
					dcomplex v1, v2, v3, v4, res1, res2, res;
					int c1, c2, c3, c4;
					c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
					c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
					c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
					c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];

					if ((c1 + c3) * (c2 + c4) > 0)
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

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;

		}
	}

	m->a_mat = new crsMatrix(*result_a_matrix);

	delete result_a_matrix;
}

void sort_matrix(Tensor_Coordinates * matrix)
{
	unsigned int key;
	unsigned int pos;
	MKL_Complex16 data_tmp;
	unsigned int int_tmp;

	for (unsigned int i = 0; i < matrix[0].k; i++)
	{
		key = matrix[0].hash[i];
		pos = i;
		for (unsigned int j = i + 1; j < matrix[0].k; j++)
		{
			if (matrix[0].hash[j] < key)
			{
				key = matrix[0].hash[j];
				pos = j;
			}
		}
		if (pos != i)
		{
			data_tmp = matrix[0].data[i];
			matrix[0].data[i] = matrix[0].data[pos];
			matrix[0].data[pos] = data_tmp;

			int_tmp = matrix[0].coord1[i];
			matrix[0].coord1[i] = matrix[0].coord1[pos];
			matrix[0].coord1[pos] = int_tmp;

			int_tmp = matrix[0].coord2[i];
			matrix[0].coord2[i] = matrix[0].coord2[pos];
			matrix[0].coord2[pos] = int_tmp;

			int_tmp = matrix[0].coord3[i];
			matrix[0].coord3[i] = matrix[0].coord3[pos];
			matrix[0].coord3[pos] = int_tmp;

			int_tmp = matrix[0].hash[i];
			matrix[0].hash[i] = matrix[0].hash[pos];
			matrix[0].hash[pos] = int_tmp;
		}
	}
}
void allocMemMat(Tensor_Coordinates_1 * mat)
{
	mat->coord2 = new unsigned  int[mat->N];
	mat->coord3 = new unsigned  int[mat->N];
	mat->data = new MKL_Complex16[mat->N];
}
void swap_row(Tensor_Coordinates_1 &mat, int i, int j)
{
	unsigned int t;
	MKL_Complex16 v;
	t = mat.coord2[i]; mat.coord2[i] = mat.coord2[j]; mat.coord2[j] = t;
	t = mat.coord3[i]; mat.coord3[i] = mat.coord3[j]; mat.coord3[j] = t;
	v = mat.data[i]; mat.data[i] = mat.data[j]; mat.data[j] = v;
}
void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat)
{
	int i = 0;
	int NN = matrix->k;
	for (i = 0; i < Nmat; i++)
	{
		mat_res[i].coord1 = i;
		mat_res[i].N = 0;
	}

	for (i = 0; i < NN; i++)
	{
		int jj = matrix->coord1[i];
		if (jj < Nmat)
		{
			mat_res[jj].N++;
		}
		else
		{
			printf("*");
		}
	}

	for (i = 0; i < Nmat; i++)
	{
		allocMemMat(mat_res + i);
		mat_res[i].N = 0;
	}

	int jj, ii, k;
	for (i = 0; i < NN; i++)
	{
		jj = matrix->coord1[i];
		ii = mat_res[jj].N;
		mat_res[jj].coord2[ii] = matrix->coord2[i];
		mat_res[jj].coord3[ii] = matrix->coord3[i];
		mat_res[jj].data[ii] = matrix->data[i];
		mat_res[jj].N++;
	}

	for (i = 0; i < Nmat; i++)
	{
		Tensor_Coordinates_1 tmp = mat_res[i];
		for (ii = 0; ii < tmp.N - 1; ii++)
		{
			for (jj = 0; jj < tmp.N - 1; jj++)
			{
				if (tmp.coord2[jj] > tmp.coord2[jj + 1])
				{
					swap_row(tmp, jj, jj + 1);
				}
				else
				{
					if (tmp.coord2[jj] == tmp.coord2[jj + 1])
					{
						if (tmp.coord3[jj] > tmp.coord3[jj + 1])
						{
							swap_row(tmp, jj, jj + 1);
						}
					}
				}
			}
		}
	}
}
void fijk_coord(Tensor_Coordinates * f_ijk, int N)
{
	unsigned int size = 5 * N * N * N - 9 * N * N - 2 * N + 6;
	f_ijk[0].data = new MKL_Complex16[size];
	memset(f_ijk[0].data, 0, size * sizeof(MKL_Complex16));
	f_ijk[0].coord1 = new unsigned int[size];
	memset(f_ijk[0].coord1, 0, size * sizeof(unsigned int));
	f_ijk[0].coord2 = new unsigned int[size];
	memset(f_ijk[0].coord2, 0, size * sizeof(unsigned int));
	f_ijk[0].coord3 = new unsigned int[size];
	memset(f_ijk[0].coord3, 0, size * sizeof(unsigned int));
	f_ijk[0].hash = new unsigned int[size];
	memset(f_ijk[0].hash, 0, size * sizeof(unsigned int));
	f_ijk[0].k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk[0].data[f_ijk[0].k].real = (i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(i), IndJ(i, j), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(i), IndS(i, j), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(i), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(i), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(i);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(i));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(i);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(i));
			f_ijk[0].k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndD(m), IndJ(i, j), IndS(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndD(m), IndS(i, j), IndJ(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(m), IndS(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(m), IndJ(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndD(m);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(m));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndD(m);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(m));
				f_ijk[0].k++;
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk[0].data[f_ijk[0].k].real = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(j), IndS(i, j), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(j), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(j));
			f_ijk[0].k++;
		}
	}

	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if (k > j)
				{
					if (k > i)
					{
						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Jjk*[Sij,Sik]  //and symmetric Ski, Jkj
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(j, k), IndS(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Jik*[Sij,Sjk]  //and symmetric Skj, Jki
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndS(j, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Sjk*[Sij,Jik]  //and symmetric Skj, Jki
						f_ijk[0].coord1[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(j, k), IndS(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Sik*[Sij,Jjk]  //and symmetric Ski, Jkj
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndJ(j, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Sjk*[Jij,Sik]  //and symmetric Ski, Skj
						f_ijk[0].coord1[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(j, k), IndJ(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Sik*[Jij,Sjk]  //and symmetric Skj, Ski
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndS(j, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Jjk*[Jij,Jik]  //and symmetric Jkj, Jki
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(j, k), IndJ(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Jik*[Jij,Jjk]  //and symmetric Jki, Jkj
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndJ(j, k));
						f_ijk[0].k++;
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndJ(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndJ(k, j));
						f_ijk[0].k++;
					}
					if (k < i)
					{
						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndS(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, i), IndS(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndJ(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, i), IndS(i, j), IndJ(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndS(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, i), IndJ(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndJ(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, i), IndJ(i, j), IndJ(k, j));
						f_ijk[0].k++;
					}
				}
			}
		}
	}
}
void dijk_coord(Tensor_Coordinates * d_ijk, int N)
{
	unsigned int size = 6 * N * N * N - (N * (21 * N + 7)) / 2 + 1;
	d_ijk[0].data = new MKL_Complex16[size];
	memset(d_ijk[0].data, 0, size * sizeof(MKL_Complex16));
	d_ijk[0].coord1 = new unsigned int[size];
	memset(d_ijk[0].coord1, 0, size * sizeof(unsigned int));
	d_ijk[0].coord2 = new unsigned int[size];
	memset(d_ijk[0].coord2, 0, size * sizeof(unsigned int));
	d_ijk[0].coord3 = new unsigned int[size];
	memset(d_ijk[0].coord3, 0, size * sizeof(unsigned int));
	d_ijk[0].hash = new unsigned int[size];
	memset(d_ijk[0].hash, 0, size * sizeof(unsigned int));
	d_ijk[0].k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndS(i, j), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(i), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndJ(i, j), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(i), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(i));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(i));
			d_ijk[0].k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(m), IndS(i, j), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(m), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(m), IndJ(i, j), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(m), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(m);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(m));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(m);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(m));
				d_ijk[0].k++;
			}
		}
	}

	for (int j = 2; j < N; j++)
	{
		int i = 0;
		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndS(i, j), IndS(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(j), IndS(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
		d_ijk[0].k++;
	}

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndS(i, j), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(j), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
			d_ijk[0].k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int z = j + 1; z < N; z++)
			{
				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(z), IndS(i, j), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(z), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(z), IndJ(i, j), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(z), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(z);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(z));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(z);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(z));
				d_ijk[0].k++;
			}
		}
	}

	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if (k > j)
				{
					if (k > i)
					{
						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Sjk*{Sij,Sik}  //and symmetric Ski, Skj
						d_ijk[0].coord1[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(j, k), IndS(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Sik*{Sij,Sjk}  //and symmetric Skj, Ski
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndS(j, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Jjk*{Sij,Jik}  //and symmetric Jkj, Jki
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(j, k), IndS(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Jik*{Sij,Jjk}  //and symmetric Jki, Jkj
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndJ(j, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);  //Jjk*{Jij,Sik}  //and symmetric Ski, Jkj
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(j, k), IndJ(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Jik*{Jij,Sjk}  //and symmetric Skj, Jki
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndS(j, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Sjk*{Jij,Jik}  //and symmetric Skj, Jki
						d_ijk[0].coord1[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(j, k), IndJ(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);  //Sik*{Jij,Jjk}  //and symmetric Ski, Jkj
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndJ(j, k));
						d_ijk[0].k++;
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndJ(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndJ(k, j));
						d_ijk[0].k++;
					}
					if (k < i)
					{
						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndS(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, i), IndS(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndJ(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, i), IndS(i, j), IndJ(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndS(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, i), IndJ(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndJ(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, i), IndJ(i, j), IndJ(k, j));
						d_ijk[0].k++;
					}
				}
			}
		}
	}


	for (int j = 2; j < N; j++)
	{
		int i = 1;
		d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndD(i), IndD(i));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(j), IndD(i));
		d_ijk[0].k++;
	}

	for (int i = 2; i < N; i++)
	{
		d_ijk[0].data[d_ijk[0].k].real = 2.0 * (1 - i) / sqrt((double)(i) * (i + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(i));
		d_ijk[0].k++;

		for (int j = i + 1; j < N; j++)
		{
			d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndD(i), IndD(i));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(j), IndD(i));
			d_ijk[0].k++;
		}
	}
}
void print_matrix(Tensor_Coordinates * matrix, FILE * f)
{
	for (unsigned int i = 0; i < matrix[0].k; i++)
	{
		fprintf(f, "%2.3i %2.3i %2.3i %lf\n", matrix[0].coord1[i] + 1, matrix[0].coord2[i] + 1, matrix[0].coord3[i] + 1, matrix[0].data[i].real);
	}
}
void free_matrix(Tensor_Coordinates * matrix)
{
	delete(matrix[0].data);
	delete(matrix[0].coord1);
	delete(matrix[0].coord2);
	delete(matrix[0].coord3);
	delete(matrix[0].hash);
	matrix[0].k = 0;
}

void init_f_d(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;

	FMatrixs * Fs = m->Fs;

	crsMatrix * Fjk, *Fkj, *Fk, *Fsum, *Fsub, *Fres;

	dcomplex add, sub;
	add.re = 1.0; add.im = 0.0;
	sub.re = -1.0; sub.im = 0.0;

	vector<pair<int, dcomplex> > * f_std;
	vector<pair<int, dcomplex> > * d_std;

	m->d_mat = new crsMatrix *[N_mat];
	m->f_mat = new crsMatrix *[N_mat];

	int i, j, k, cnt;

	Fjk = new crsMatrix;
	Fkj = new crsMatrix;

	Fsum = new crsMatrix;
	Fsub = new crsMatrix;

	Fres = new crsMatrix;

	for (i = 0; i < N_mat; i++)
	{
		f_std = new vector<pair<int, dcomplex> >[N_mat];
		d_std = new vector<pair<int, dcomplex> >[N_mat];

		for (j = 0; j < N_mat; j++)
		{
			for (k = 0; k < N_mat; k++)
			{

				Fk = new crsMatrix(*(Fs->F[k + 1]));

				SparseMKLMult(*(Fs->F[j + 1]), *Fk, *Fjk, true);
				SparseMKLMult(*Fk, *(Fs->F[j + 1]), *Fkj, true);
				delete Fk;

				SparseMKLAdd(*Fjk, add, *Fkj, *Fsum, true);
				SparseMKLAdd(*Fjk, sub, *Fkj, *Fsub, true);


				SparseMKLMult(*(Fs->F[i + 1]), *Fsub, *Fres, true);
				cnt = trace_struct(*Fres);
				if (cnt > 0)
				{
					dcomplex val, tr;
					tr = trace(*Fres);
					val.re = tr.im;
					val.im = -tr.re;
					f_std[j].push_back(make_pair(k, val));
				}

				SparseMKLMult(*(Fs->F[i + 1]), *Fsum, *Fres, true);
				cnt = trace_struct(*Fres);
				if (cnt > 0)
				{
					dcomplex val;
					val = trace(*Fres);
					d_std[j].push_back(make_pair(k, val));
				}
			}
		}

		m->f_mat[i] = stdToCrs(f_std, N_mat);
		m->d_mat[i] = stdToCrs(d_std, N_mat);

		//    printMatrixVal(m->f_mat[i]);
		//    printMatrixVal(m->d_mat[i]);

		delete[] f_std;
		delete[] d_std;
	}
	delete Fjk;
	delete Fkj;

	delete Fsum;
	delete Fsub;

	delete Fres;

}
crsMatrix * TensorToCrs(Tensor_Coordinates &t_ijk, int s, int f, int N)
{
	crsMatrix *  mat = new crsMatrix(N, f - s);
	int i = 0;
	dcomplex val;

	for (i = 0; i < f - s; i++)
	{
		val.re = t_ijk.data[s + i].real;
		val.im = t_ijk.data[s + i].imag;
		mat->Value[i] = val;
		mat->Col[i] = t_ijk.coord3[s + i];
	}

	for (i = 0; i < N + 1; i++)
	{
		mat->RowIndex[i] = 0;
	}

	for (i = 0; i < f - s; i++)
	{
		mat->RowIndex[t_ijk.coord2[s + i] + 1]++;
	}

	int cnt = 0;

	for (i = 0; i < N; i++)
	{
		mat->RowIndex[i] = cnt;
		cnt += mat->RowIndex[i + 1];
	}
	mat->RowIndex[i] = cnt;

	return mat;
}
crsMatrix * TensorToCrs(Tensor_Coordinates_1 &t_ijk, int s, int f, int N)
{
	crsMatrix *  mat = new crsMatrix(N, f - s);
	int i = 0;
	dcomplex val;

	for (i = 0; i < f - s; i++)
	{
		val.re = t_ijk.data[s + i].real;
		val.im = t_ijk.data[s + i].imag;
		mat->Value[i] = val;
		mat->Col[i] = t_ijk.coord3[s + i];
	}

	for (i = 0; i < N + 1; i++)
	{
		mat->RowIndex[i] = 0;
	}

	for (i = 0; i < f - s; i++)
	{
		mat->RowIndex[t_ijk.coord2[s + i] + 1]++;
	}

	int cnt = 0;

	for (i = 0; i < N; i++)
	{
		mat->RowIndex[i] = cnt;
		cnt += mat->RowIndex[i + 1];
	}
	mat->RowIndex[i] = cnt;

	return mat;
}
void init_f_d_valentin(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	int i;

	m->d_mat = new crsMatrix *[N_mat];
	m->f_mat = new crsMatrix *[N_mat];

	//  vector<pair<int , dcomplex> > * f_std;
	//  vector<pair<int , dcomplex> > * d_std;

	Tensor_Coordinates f_ijk;
	Tensor_Coordinates d_ijk;
	Tensor_Coordinates_1 * f_1 = new Tensor_Coordinates_1[N_mat];
	Tensor_Coordinates_1 * d_1 = new Tensor_Coordinates_1[N_mat];

	fijk_coord(&f_ijk, N + 1);
	sort_matrix(&f_ijk, f_1, N_mat);
	//  sort_matrix(&f_ijk);

	dijk_coord(&d_ijk, N + 1);
	sort_matrix(&d_ijk, d_1, N_mat);
	//  sort_matrix(&d_ijk);

	int k, s, f;

	//  k = 0;
	for (i = 0; i < N_mat; i++)
	{
		//    s = k;
		//    while(f_ijk.coord1[k] == i)
		//    {
		//      k++;
		//    }
		//    f = k;
		//    m->f_mat[i] = TensorToCrs(f_ijk, s, f, N_mat);
		m->f_mat[i] = TensorToCrs(f_1[i], 0, f_1[i].N, N_mat);
	}

	//  k = 0;
	for (i = 0; i < N_mat; i++)
	{
		//    s = k;
		//    while(d_ijk.coord1[k] == i)
		//    {
		//      k++;
		//    }
		//    f = k;
		//    m->d_mat[i] = TensorToCrs(d_ijk, s, f, N_mat);
		m->d_mat[i] = TensorToCrs(d_1[i], 0, d_1[i].N, N_mat);
	}

	for (i = 0; i < N_mat; i++)
	{
		delete[] f_1[i].coord2;
		delete[] f_1[i].coord3;
		delete[] f_1[i].data;
		delete[] d_1[i].coord2;
		delete[] d_1[i].coord3;
		delete[] d_1[i].data;
	}
	delete[] f_1;
	delete[] d_1;

	free_matrix(&f_ijk);
	free_matrix(&d_ijk);
}

void transpFs(Model *m)
{
	int N_mat = m->N_mat;
	crsMatrix **f_mat = m->f_mat;
	crsMatrix **f_H_mat = new crsMatrix*[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		f_H_mat[i] = new crsMatrix(*(f_mat[i]));
		Transpose(*(f_mat[i]), *(f_H_mat[i]), false);
		//    printMatrixVal(f_H_mat[i]);
	}
	m->f_H_mat = f_H_mat;
}

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
	for (int i = 1; i <= N_mat; i++)
	{
		QEs->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((hE[i].re != 0.0) || (hE[i].im != 0.0))
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
	for (int i = 1; i <= N_mat; i++)
	{
		QUs->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((hU[i].re != 0.0) || (hU[i].im != 0.0))
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
	for (int i = 1; i <= N_mat; i++)
	{
		QJs->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((hJ[i].re != 0.0) || (hJ[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*QJs, hJ[i], *(f_mat[i]), *resSum);
			delete QJs;
			QJs = resSum;
		}
	}

	m->QJs = QJs;
}

void calcKs(Model *m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;
	int N_mat = m->N_mat;
	crsMatrix *Ks_tmp;
	crsMatrix *FsT;
	dcomplex  *Ks = m->Ks;
	crsMatrix *As = m->a_mat;
	crsMatrix **Fs = m->f_mat;
	crsMatrix *AsT;

	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re = 0.0;
		Ks[i].im = 0.0;
	}

	//  printMatrixVal(As);
	for (int i = 0; i < N_mat; i++)
	{
		AsT = As;
		//AsT = new crsMatrix(*(As));
		//Transpose(*(As), *AsT);
		FsT = Fs[i];
		FsT = new crsMatrix(*(Fs[i]));
		Transpose(*(Fs[i]), *FsT, false);
		//    printMatrixVal(FsT);
		for (int j = 0; j < N_mat; j++)
		{
			int ii, jj, m1, m2;
			ii = As->RowIndex[i];
			m1 = As->RowIndex[i + 1];
			jj = FsT->RowIndex[j];
			m2 = FsT->RowIndex[j + 1];

			while ((ii < m1) && (jj < m2))
			{
				if (AsT->Col[ii] < FsT->Col[jj])
				{
					ii++;
				}
				else
				{
					if (AsT->Col[ii] > FsT->Col[jj])
					{
						jj++;
					}
					else
					{
						dcomplex as, fs;
						as = AsT->Value[ii];
						fs = FsT->Value[jj];
						Ks[j].re += as.re * fs.re - as.im * fs.im;
						Ks[j].im += as.re * fs.im + as.im * fs.re;

						ii++;
						jj++;
					}
				}
			}

		}
		delete FsT;
		//delete AsT;
	}

	dcomplex val;
	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re *= (-cp.g) / (N + 1);
		Ks[i].im *= (cp.g) / (N + 1);
		val = Ks[i];
		Ks[i].re = val.im;
		Ks[i].im = val.re;
	}
}

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
	for (int i = 1; i <= N_mat; i++)
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

	for (int i = 0; i < N_mat; i++)
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
	for (int i = 0; i < N_mat; i++)
	{
		toZeroBase(*(f_mat[i]));
		toZeroBase(*(d_mat[i]));
		toZeroBase(*(f_H_mat[i]));
	}
	//  delete FkT;
	//  delete FiT;

	double mm = -cp.g / 4.0;
	for (int i = 0; i < Rs->NZ; i++)
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
	crsMatrix * MatSDi, *MatSDk;
	crsMatrix * M1, *M2, *Rsum;

	crsMatrix * Rs_tmp;
	crsMatrix * Rs = new crsMatrix(m->N_mat, 1);

	Rs->Col[0] = 1;
	Rs->RowIndex[0] = 1;
	for (int i = 1; i <= N_mat; i++)
	{
		Rs->RowIndex[i] = 1;
	}

	crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_H_mat = m->f_H_mat;
	crsMatrix ** d_mat = m->d_mat;
	crsMatrix *  As = m->a_mat;

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

	for (int i = start; i < finish; i++)
	{
		SparseMKLAddOne(*(f_mat[i]), sum_i, *(d_mat[i]), *MatSDi, true);

		int st = As->RowIndex[i];
		int fn = As->RowIndex[i + 1];

		for (int j = st; j < fn; j++)
		{
			int k = As->Col[j];
			SparseMKLAddOne(*(f_mat[k]), sum_i, *(d_mat[k]), *MatSDk, true);

			for (int conj_i = 0; conj_i < MatSDk->NZ; conj_i++)
				MatSDk->Value[conj_i].im = -MatSDk->Value[conj_i].im;

			SparseMKLMultOne(*(f_H_mat[k]), *MatSDi, *M1, true);
			SparseMKLMultOne(*(f_H_mat[i]), *MatSDk, *M2, true);

			SparseMKLAddOne(*M1, sum, *M2, *Rsum, true);

			dcomplex val1 = As->Value[j];
			for (int ii = 0; ii < Rsum->NZ; ii++)
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

void calcGs(Model *m)
{
	int N_mat = m->N_mat;
	crsMatrix * Gs = new crsMatrix();

	crsMatrix * subSum1 = new crsMatrix();
	crsMatrix * subSum2 = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->QEs), sum, *(m->QUs), *subSum1);
	SparseMKLAdd(*(m->QJs), sum, *(m->Rs), *subSum2);

	SparseMKLAdd(*(subSum1), sum, *(subSum2), *Gs);
	m->Gs = Gs;
}

int allocComplexMatrix(int n, double ** re, double ** im)
{
	double * mem;
	mem = (double *)malloc(sizeof(double) * n * n);
	if (mem == NULL)
	{
		return -1;
	}
	(*re) = mem;

	mem = (double *)malloc(sizeof(double) * n * n);
	if (mem == NULL)
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
int genNormalDistributedElemets(int n, double * re, double * im)
{
	return genNormalDistributedElemets(n, n, re, im);
}
int genNormalDistributedElemets(int n1, int n2, double * re, double * im)
{
	VSLStreamStatePtr stream;
	int i, errcode;
	double a = 0.0, sigma = 1.0;

	/***** Initialize *****/
	errcode = vslNewStream(&stream, BRNG, SEED);
	if (errcode != VSL_STATUS_OK) return errcode;

	errcode = vdRngGaussian(METHOD, stream, n1 * n2, re, a, sigma);
	if (errcode != VSL_STATUS_OK) return errcode;

	errcode = vdRngGaussian(METHOD, stream, n1 * n2, im, a, sigma);
	if (errcode != VSL_STATUS_OK) return errcode;

	/***** Deinitialize *****/
	errcode = vslDeleteStream(&stream);
	if (errcode != VSL_STATUS_OK) return errcode;

	return 0;
}
int createHermitianMatrix(int n, dcomplex *a)
{
	double *re1, *im1, *re2, *im2;
	int errcode, i, j;
	errcode = allocComplexMatrix(n, &re1, &im1);
	if (errcode != 0) return -3;

	errcode = allocComplexMatrix(n, &re2, &im2);
	if (errcode != 0) return -4;

	genNormalDistributedElemets(n, re1, im1);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			re2[i*n + j] = re1[j*n + i];
			im2[i*n + j] = -im1[j*n + i];
		}
	}

	double sre = 0.0, sim = 0.0;

	double two = 2.0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			re1[i*n + j] = (re1[i*n + j] + re2[i*n + j]) / two;
			im1[i*n + j] = (im1[i*n + j] + im2[i*n + j]) / two;
			sre += fabs(re1[i*n + j]);
			sim += fabs(im1[i*n + j]);
		}
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			re1[i*n + j] /= sre;
			im1[i*n + j] /= sim;
		}
	}

	freeComplexMatrix(re2, im2);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
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
	double sum = 0.0;

	errcode = allocComplexMatrix(n, &re, &im);
	if (errcode != 0) return -3;

	genNormalDistributedElemets(n, 1, re, im);
	for (i = 0; i < n; i++)
	{
		sum += re[i] * re[i] + im[i] * im[i];
	}

	for (i = 0; i < n; i++)
	{
		re[i] = sqrt(1.0 / n / 2);
		im[i] = sqrt(1.0 / n / 2);
	}
	sum = 0.0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i*n + j].re = re[i] * re[j] + im[i] * im[j];
			a[i*n + j].im = re[i] * im[j] - im[i] * re[j];
		}
		sum += a[i*n + i].re;
	}

	printf("Init Tr(Rho): %lf\n", sum);

	freeComplexMatrix(re, im);

	return 0;
}

void init_conditions(dcomplex *mtx, ConfigParam &cp, ConfigData &cd)
{
	for (int s_id_1 = 0; s_id_1 < cd.Ns; s_id_1++)
	{
		for (int s_id_2 = 0; s_id_2 < cd.Ns; s_id_2++)
		{
			mtx[s_id_1 * cd.Ns + s_id_2].re = 0.0;
			mtx[s_id_1 * cd.Ns + s_id_2].im = 0.0;
		}
	}

	if (cp.int_ist == 0)
	{
		mtx[cp.int_isi * cd.Ns + cp.int_isi].re = 1.0;
		mtx[cp.int_isi * cd.Ns + cp.int_isi].im = 0.0;
	}
	else if (cp.int_ist == 1)
	{
		for (int s_id_1 = 0; s_id_1 < cd.Ns; s_id_1++)
		{
			mtx[s_id_1 * cd.Ns + s_id_1].re = 1.0 / cd.Ns;
			mtx[s_id_1 * cd.Ns + s_id_1].im = 0.0;
		}
	}
	else
	{
		stringstream msg;
		msg << "wrong int_ist value: " << cp.int_ist << endl;
		Error(msg.str());
	}
}

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res)
{
	int i, j, s, f;
	for (i = 0; i < mat->N; i++)
	{
		s = mat->RowIndex[i];
		f = mat->RowIndex[i + 1];
		res[i].re = 0.0;
		res[i].im = 0.0;
		for (j = s; j < f; j++)
		{
			dcomplex v1 = mat->Value[j];
			dcomplex v2 = x[mat->Col[j]];
			res[i].re += v1.re * v2.re;// - v1.im * v2.im;
			res[i].im = 0.0;//+= v1.re * v2.im + v1.im * v2.re;
		}
	}
}
void calcVectValue(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp1, dcomplex * tmp2)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * Gs = m->Gs;
	crsMatrix * QEs = m->QEs;
	dcomplex  * Ks = m->Ks;
	double T = 0.0;
	double A0 = 0.0;
	double w = 0.0;

	multMatVec(Gs, x, tmp1);
	multMatVec(QEs, x, tmp2);

	for (i = 0; i < N_mat; i++)
	{
		res[i].re = (tmp1[i].re - Ks[i].re) * h;
		res[i].im = 0.0;
	}
}
dcomplex calcDiffIter(Model *m)
{
	dcomplex max_diff, diff;
	max_diff.re = 0.0;
	max_diff.im = 0.0;
	for (int i = 0; i < m->N_mat; i++)
	{
		diff.re = fabs(m->prevRhoF[i].re - m->RhoF[i].re);
		diff.im = fabs(m->prevRhoF[i].im - m->RhoF[i].im);
		if (max_diff.re < diff.re)max_diff.re = diff.re;
		if (max_diff.im < diff.im)max_diff.im = diff.im;
	}

	return max_diff;
}
void calcODE(Model *m, IntData &int_data, ConfigData &cd, ConfigParam &cp)
{
	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = int_data.k1;
	dcomplex * k2 = int_data.k2;
	dcomplex * k3 = int_data.k3;
	dcomplex * k4 = int_data.k4;
	dcomplex * val = int_data.val;
	dcomplex * tmp1 = int_data.tmp1;
	dcomplex * tmp2 = int_data.tmp2;
	dcomplex * tmp3 = int_data.tmp3;


	double curr_time = 0.0;
	cout << endl << "dump_time: " << curr_time << endl;
	clearRho(m);
	calcRho(m);
	characteristics(m, cd, cp, false);


	int curr_dump = 1;
	double curr_dump_time = int_data.dump_times[curr_dump];
	double curr_step = cp.int_h;

	while (curr_dump < cp.int_dn + 2)
	{
		while (fabs(curr_time + curr_step - curr_dump_time) > eps &&  curr_time + curr_step < curr_dump_time)
		{
			calcVectValue(curr_step, m, RhoF, k1, tmp1, tmp2);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue(curr_step, m, val, k2, tmp1, tmp2);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue(curr_step, m, val, k3, tmp1, tmp2);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue(curr_step, m, val, k4, tmp1, tmp2);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}

			curr_time += curr_step;
		}

		// Last iteration and dumping
		curr_step = curr_dump_time - curr_time;

		calcVectValue(curr_step, m, RhoF, k1, tmp1, tmp2);
		for (i = 0; i < N_mat; i++)
		{
			val[i].re = RhoF[i].re + k1[i].re / 2.0;
			val[i].im = RhoF[i].im + k1[i].im / 2.0;
		}
		calcVectValue(curr_step, m, val, k2, tmp1, tmp2);
		for (i = 0; i < N_mat; i++)
		{
			val[i].re = RhoF[i].re + k2[i].re / 2.0;
			val[i].im = RhoF[i].im + k2[i].im / 2.0;
		}
		calcVectValue(curr_step, m, val, k3, tmp1, tmp2);
		for (i = 0; i < N_mat; i++)
		{
			val[i].re = RhoF[i].re + k3[i].re;
			val[i].im = RhoF[i].im + k3[i].im;
		}
		calcVectValue(curr_step, m, val, k4, tmp1, tmp2);

		for (i = 0; i < N_mat; i++)
		{
			RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
			RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
		}

		curr_time = int_data.dump_times[curr_dump];

		cout << endl << "dump_time: " << curr_time << endl;
		clearRho(m);
		calcRho(m);
		characteristics(m, cd, cp, true);

		curr_dump++;
		curr_dump_time = int_data.dump_times[curr_dump];
		curr_step = cp.int_h;
	}
}

dcomplex multMatCRS_tr(dcomplex *a, crsMatrix *b)
{
	int N = b->N;

	dcomplex c;
	c.re = 0.0;
	c.im = 0.0;

	crsMatrix *tB = new crsMatrix(*b);
	Transpose(*b, *tB, false);

	int i, j, jj, k, s, f;
	for (i = 0; i < N; i++)
	{
		s = tB->RowIndex[i];
		f = tB->RowIndex[i + 1];
		for (k = s; k < f; k++)
		{
			j = tB->Col[k];
			c.re += a[i * N + j].re * tB->Value[k].re - a[i * N + j].im * tB->Value[k].im;
			c.im += a[i * N + j].re * tB->Value[k].im + a[i * N + j].im * tB->Value[k].re;
		}
	}

	delete tB;

	return c;
}
void initRhoODE(Model *m, ConfigData &cd, ConfigParam &cp)
{
	int N = m->N;
	int N_mat = m->N_mat;
	FMatrixs  * Fs = m->Fs;
	dcomplex  * RhoF = m->RhoF;
	dcomplex * psi = new dcomplex[(N + 1)*(N + 1)];
	init_conditions(psi, cp, cd);

	for (int i = 0; i < N_mat; i++)
	{
		RhoF[i] = multMatCRS_tr(psi, Fs->F[i + 1]);
	}

	delete[] psi;
}

void linSolv(Model *m)
{
	crsMatrix * Gs = m->Gs;
	dcomplex  * Ks = m->Ks;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;

	for (int i = 0; i < Gs->NZ; i++)
		Gs->Col[i]++;
	for (int j = 0; j <= N_mat; j++)
	{
		Gs->RowIndex[j]++;
	}


	MKL_INT mtype = 13; /* Real unsymmetric matrix */

	MKL_INT nrhs = 1; /* Number of right hand sides. */
					  /* Internal solver memory pointer pt, */
					  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
					  /* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum; /* Double dummy */
	MKL_INT idum; /* Integer dummy. */
				  /* -------------------------------------------------------------------- */
				  /* .. Setup Pardiso control parameters. */
				  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
				  /* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
			   /* -------------------------------------------------------------------- */
			   /* .. Initialize the internal solver memory pointer. This is only */
			   /* necessary for the FIRST call of the PARDISO solver. */
			   /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, Gs->Value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	//  printf("\nReordering completed ... ");
	//  printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	//  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, Gs->Value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	//  printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
				  /* Set right hand side to one. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, Gs->Value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, Ks, RhoF, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, &ddum, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	for (int i = 0; i < Gs->NZ; i++)
		Gs->Col[i]--;
	for (int j = 0; j <= N_mat; j++)
	{
		Gs->RowIndex[j]--;
	}
}
void linSolvCheck(Model *m)
{
	crsMatrix * Gs = m->Gs;
	dcomplex  * Ks = m->Ks;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;

	dcomplex * Ks_r = new dcomplex[N_mat];

	multMatVec(Gs, RhoF, Ks_r);

	double m1, m2;

	m1 = fabs(Ks[0].re - Ks_r[0].re);
	m2 = fabs(Ks[0].im - Ks_r[0].im);

	for (int i = 1; i < N_mat; i++)
	{
		if (m1 < fabs(Ks[0].re - Ks_r[0].re)) m1 = fabs(Ks[0].re - Ks_r[0].re);
		if (m2 < fabs(Ks[0].im - Ks_r[0].im)) m2 = fabs(Ks[0].im - Ks_r[0].im);
	}

	printf("diff (%1.16lf + i %1.16lf) \n", m1, m2);

	delete[] Ks_r;
}
void linSolvReal(Model *m)
{
	crsMatrix * Gs = m->Gs;
	dcomplex  * Ks = m->Ks;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;
	double * value, *res;
	value = new double[Gs->NZ];
	res = new double[Gs->N];

	for (int i = 0; i < Gs->NZ; i++)
	{
		Gs->Col[i]++;
		value[i] = Gs->Value[i].re;
	}
	for (int j = 0; j <= N_mat; j++)
	{
		Gs->RowIndex[j]++;
	}


	MKL_INT mtype = 11; /* Real unsymmetric matrix */

	MKL_INT nrhs = 1; /* Number of right hand sides. */
					  /* Internal solver memory pointer pt, */
					  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
					  /* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum; /* Double dummy */
	MKL_INT idum; /* Integer dummy. */
				  /* -------------------------------------------------------------------- */
				  /* .. Setup Pardiso control parameters. */
				  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
				  /* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
			   /* -------------------------------------------------------------------- */
			   /* .. Initialize the internal solver memory pointer. This is only */
			   /* necessary for the FIRST call of the PARDISO solver. */
			   /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	//  printf("\nReordering completed ... ");
	//  printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	//  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	//  printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
				  /* Set right hand side to one. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, value, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, Ks, res, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	for (int i = 0; i < Gs->N; i++)
	{
		RhoF[i].re = res[i];
		RhoF[i].im = 0.0;
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N_mat, &ddum, Gs->RowIndex, Gs->Col, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	for (int i = 0; i < Gs->NZ; i++)
		Gs->Col[i]--;
	for (int j = 0; j <= N_mat; j++)
	{
		Gs->RowIndex[j]--;
	}
	delete[]res;
	delete[]value;
}

void calcRho(Model *m)
{
	FMatrixs  * Fs = m->Fs;
	crsMatrix * Rho = new crsMatrix(*(Fs->F[0]));
	crsMatrix * Rho_tmp = new crsMatrix(*(Fs->F[0]));
	crsMatrix * tmp;
	dcomplex  * RhoF = m->RhoF;
	int N_mat = m->N_mat;
	int N = m->N;

	int i;
	for (i = 0; i < Rho->NZ; i++)
	{
		Rho->Value[i].re /= N + 1;
		Rho->Value[i].im /= N + 1;
	}

	for (int i = 0; i < N_mat; i++)
	{
		SparseMKLAdd(*Rho, RhoF[i], *(Fs->F[i + 1]), *Rho_tmp, true);
		tmp = Rho; Rho = Rho_tmp; Rho_tmp = tmp;
	}
	delete Rho_tmp;

	m->Rho = Rho;
}

void clearRho(Model *m)
{
	if (m->Rho)
	{
		delete m->Rho;
		m->Rho = NULL;
	}
}


void characteristics(Model *m, ConfigData &cd, ConfigParam &cp, bool append)
{
	double * characteristics = new double[1];

	MKL_Complex16 * rho_in_d = new MKL_Complex16[cd.Ns * cd.Ns];

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
		{
			rho_in_d[state_id_1 * cd.Ns + state_id_2].real = 0.0;
			rho_in_d[state_id_1 * cd.Ns + state_id_2].imag = 0.0;
		}
	}

	for (int i = 0; i < m->Rho->N; i++)
	{
		for (int k = m->Rho->RowIndex[i]; k < m->Rho->RowIndex[i + 1]; k++)
		{
			rho_in_d[i * cd.Ns + m->Rho->Col[k]].real = m->Rho->Value[k].re;
			rho_in_d[i * cd.Ns + m->Rho->Col[k]].imag = m->Rho->Value[k].im;
		}
	}

	// ######## imbalance ########
	int Nss = pow(2, cp.Nc);

	MKL_Complex16 * rho_xtd = new MKL_Complex16[Nss * Nss];

	for (int xtd_st_id_1 = 0; xtd_st_id_1 < Nss; xtd_st_id_1++)
	{
		for (int xtd_st_id_2 = 0; xtd_st_id_2 < Nss; xtd_st_id_2++)
		{
			rho_xtd[xtd_st_id_1 * Nss + xtd_st_id_2].real = 0.0;
			rho_xtd[xtd_st_id_1 * Nss + xtd_st_id_2].imag = 0.0;
		}
	}

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
		{
			int xtd_st_id_1 = cd.id_to_x[state_id_1];
			int xtd_st_id_2 = cd.id_to_x[state_id_2];

			rho_xtd[xtd_st_id_1 * Nss + xtd_st_id_2].real = rho_in_d[state_id_1 * cd.Ns + state_id_2].real;
			rho_xtd[xtd_st_id_1 * Nss + xtd_st_id_2].imag = rho_in_d[state_id_1 * cd.Ns + state_id_2].imag;
		}
	}

	int red_size = pow(2, cp.Nc / 2);

	MKL_Complex16 * rho_red = new MKL_Complex16[red_size * red_size];

	for (int red_st_id_1 = 0; red_st_id_1 < red_size; red_st_id_1++)
	{
		for (int red_st_id_2 = 0; red_st_id_2 < red_size; red_st_id_2++)
		{
			MKL_Complex16 sum = { 0.0, 0.0 };

			for (int sum_id = 0; sum_id < red_size; sum_id++)
			{
				int rho_id_1 = red_st_id_1 * red_size + sum_id;
				int rho_id_2 = red_st_id_2 * red_size + sum_id;

				sum.real += rho_xtd[rho_id_1 * Nss + rho_id_2].real;
				sum.imag += rho_xtd[rho_id_1 * Nss + rho_id_2].imag;
			}

			rho_red[red_st_id_1 * red_size + red_st_id_2].real = sum.real;
			rho_red[red_st_id_1 * red_size + red_st_id_2].imag = sum.imag;
		}
	}

	double eps_eval = 1.0e-8;

	double * red_evals = new double[red_size];
	for (int red_st_id = 0; red_st_id < red_size; red_st_id++)
	{
		red_evals[red_st_id] = 0.0;
	}

	int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'N', 'L', red_size, rho_red, red_size, red_evals);

	double ee = 0.0; // entanglement entropy

	for (int red_st_id = 0; red_st_id < red_size; red_st_id++)
	{
		if (fabs(red_evals[red_st_id]) > eps_eval)
		{
			ee -= red_evals[red_st_id] * log(red_evals[red_st_id]);
		}	
	}

	characteristics[0] = ee;

	string characteristics_fn = cp.path + "entropy" + file_name_suffix(cp, 4);
	write_double_data(characteristics_fn, characteristics, 1, 16, append);

	delete[] red_evals;
	delete[] rho_red;
	delete[] rho_xtd;

	// ######## imbalance ########
	double * n_part = new double[cd.Nc];

	for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
	{
		n_part[cell_id] = 0.0;
	}

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		vector<int> vb = convert_int_to_vector_of_bits(cd.id_to_x[state_id], cd.Nc);

		for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
		{
			n_part[cell_id] += rho_in_d[state_id * cd.Ns + state_id].real * double(vb[cell_id]);
		}
	}

	double sum_odd = 0.0;
	double sum_even = 0.0;
	double sum_all = 0.0;


	for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
	{
		if (cell_id % 2 == 0)
		{
			sum_odd += n_part[cell_id];
		}
		else
		{
			sum_even += n_part[cell_id];
		}
	}
	sum_all = sum_even + sum_odd;

	double imbalance = (sum_odd - sum_even) / sum_all;

	delete[] n_part;

	characteristics[0] = imbalance;

	characteristics_fn = cp.path + "imbalance" + file_name_suffix(cp, 4);
	write_double_data(characteristics_fn, characteristics, 1, 16, append);

	delete[] characteristics;
	delete[] rho_in_d;
}