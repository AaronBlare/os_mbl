#include "testODE.h"
#include "calcODE.h"
#include "CalcRho.h"
#include <string.h>
#include <math.h>
#include "Matrix_op.h"

void testODE(Model *m)
{
	int N_mat = m->N_mat;
	int N = m->N;
	//int NSTEP = m->conf.NSTEP;
	int NSTEP = 100;
	double *diag = new double[m->Rho->N];

	crsMatrix * Rho_T;

	char fileName1[512];
	char fileName2[512];
	sprintf(fileName1, "diag.txt");
	sprintf(fileName2, "diag_tr.txt");


	FILE *f = fopen(fileName1, "w");
	FILE *f_tr = fopen(fileName2, "w");
	double sum;
	for (int i = 0; i < NSTEP * 5; i += 10)
	{
		//calcODE(m, m->conf.h, 10, m->conf.N_T * m->conf.T + i * m->conf.h);
		calcODE(m, 0.01, 10, 100 * 1 + i * 0.1);
		Rho_T = m->Rho;
		calcRho(m);
		AbsMatrixDiagVal(Rho_T, diag);
		sum = 0.0;
		for (int j = 0; j < Rho_T->N; j++)
		{
			sum += diag[j];
			fprintf(f, "%lf ", diag[j]);
		}
		fprintf(f, "\n");
		fprintf(f_tr, "%lf\n", sum);
		delete Rho_T;
	}
	fclose(f);
	fclose(f_tr);
	delete[]diag;
	//crsMatrix * Rho_T;
	//crsMatrix * Rho_T_div_4;
	//crsMatrix * Rho_T_div_2;
	//crsMatrix * Rho_T_div_3_4;
	//crsMatrix * Rho_2T;
	//
	//crsMatrix * div = new crsMatrix();
	//
	//Rho_T = m->Rho;
	//
	//calcODE(m, m->conf.h, m->conf.NSTEP / 4);
	//calcRho(m);
	//Rho_T_div_4 = m->Rho;
	//
	//calcODE(m, m->conf.h, m->conf.NSTEP / 4);
	//calcRho(m);
	//Rho_T_div_2 = m->Rho;
	//
	//calcODE(m, m->conf.h, m->conf.NSTEP / 4);
	//calcRho(m);
	//Rho_T_div_3_4 = m->Rho;
	//
	//calcODE(m, m->conf.h, m->conf.NSTEP / 4);
	//calcRho(m);
	//Rho_2T = m->Rho;
	//
	//dcomplex betta;
	//betta.re = -1;
	//betta.im = 0;
	//SparseMKLAdd(*Rho_T, betta, *Rho_2T, *div, true); 
	//
	//double max_re = 0.0;
	//double max_im = 0.0;
	//for(int i = 0; i < div->NZ; i++)
	//{
	//  if(max_re < abs(div->Value[i].re))max_re = abs(div->Value[i].re);
	//  if(max_im < abs(div->Value[i].re))max_re = abs(div->Value[i].re);
	//}
	//
	//printf("dif_step: %lf %lf \n", max_re, max_im);
	//
	//saveAbsMatrixVal("absRho_T_div_4.txt", Rho_T_div_4);
	//saveAbsMatrixVal("absRho_T_div_2.txt", Rho_T_div_2);
	//saveAbsMatrixVal("absRho_T_div_3_4.txt", Rho_T_div_3_4);
	//saveAbsMatrixVal("absRho_2T.txt", Rho_2T);
	//
	//delete div;
	//delete Rho_T;
	//delete Rho_T_div_2;
}
