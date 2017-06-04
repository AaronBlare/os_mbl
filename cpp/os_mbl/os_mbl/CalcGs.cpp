#include "CalcGs.h"

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
