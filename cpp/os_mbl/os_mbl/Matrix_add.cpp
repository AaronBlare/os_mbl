#include "Matrix_op.h"
#include <stdio.h>
#include "mkl.h"

// Принимает 2 квадратных матрицы в формате CRS (3 массива, индексация с нуля)
// Возвращает C = A * B, C - в формате CRS (3 массива, индексация с нуля)
//                       Память для C в начале считается не выделенной
// Возвращает признак успешности операции: 0 - ОК, 1 - не совпадают размеры (N)
// Возвращает время работы
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

// Принимает 2 квадратных матрицы в формате CRS (3 массива, индексация с нуля)
// Возвращает C = A * B, C - в формате CRS (3 массива, индексация с нуля)
//                       Память для C в начале считается не выделенной
// Возвращает признак успешности операции: 0 - ОК, 1 - не совпадают размеры (N)
// Возвращает время работы
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


// Принимает 2 квадратных матрицы в формате CRS (3 массива, индексация с нуля)
// Возвращает C = A * B, C - в формате CRS (3 массива, индексация с нуля)
//                       Память для C в начале считается не выделенной
// Возвращает признак успешности операции: 0 - ОК, 1 - не совпадают размеры (N)
// Возвращает время работы
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
