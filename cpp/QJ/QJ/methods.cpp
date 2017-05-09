#include "header.h"
#include "output.h"
#include "omp.h"
#include <iostream>
#include <mkl.h>
#include <complex>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

#define IND(a,b)  (a) * N + (b)
#define IND2(a,b) (a) * N * N + (b)


MKL_Complex16 ZERO={0,0}, ONE={1,0}, I={0,1};
double t=0.0;
MKL_Complex16 * matrix_omp;
int max_steps = 16;
int index_kT[32];
int gkT;

MKL_Complex16 * phi_global;
double * etas;
MKL_Complex16 * rho_omp;
double * abs_rho_diag_omp;

MKL_Complex16 * rho_curr;
MKL_Complex16 * rho_curr_and;
MKL_Complex16 * hamiltonian;
MKL_Complex16 * ev_hamiltonian;
MKL_Complex16 * ev_t_hamiltonian;

int * state_in_bin;
double entropy;
double imbalance;
double ipr;
double entropy_and;
double imbalance_and;
double ipr_and;
int period;


inline MKL_Complex16 Complex_mul(MKL_Complex16 a, double b)
{
	MKL_Complex16 res = {0.0, 0.0};
	res.real = a.real * b;
	res.imag = a.imag * b;
	return res;
}

inline MKL_Complex16 Complex_scalar_mul(MKL_Complex16 * a, MKL_Complex16 * b, int N)
{
	MKL_Complex16 res = {0.0, 0.0};
	for (int i = 0; i < N; i++)
	{
		res.real += a[i].real * b[i].real + a[i].imag * b[i].imag;
		res.imag += b[i].real * a[i].imag - a[i].real * b[i].imag;
	}
	return res;
}


int if_norm(MKL_Complex16 * phi, double * eta, int N)
{
	double norm=0.;
	for(int i=0; i<N; i++)
	{
		norm+=phi[i].real*phi[i].real+phi[i].imag*phi[i].imag;
	}
	if((norm)>eta[0])
		return 0;
	else
	{
		return 1;
	}
}

double norm_vector2(MKL_Complex16 * phi, int N)
{
	double norm=0.;
	for(int i=0; i<N; i++)
	{
		norm+=phi[i].real*phi[i].real+phi[i].imag*phi[i].imag;
	}
	return norm;
}

void print_complex_array(MKL_Complex16 * tmp, int N)
{
	printf("\n");
	for (int i = 0; i < N; i++)
	{	
		printf("%d: %0.2le %0.2le\n", i, tmp[i].real, tmp[i].imag);
	}
}

void print_double_array(double * tmp, int N)
{
	printf("\n");
	for (int i = 0; i < N; i++)
	{	
		printf("%d: %0.2le\n", i, tmp[i]);
	}
}

void init_data(int N,
			   int num_trajectories,
			   int num_omp_threads)
{
	phi_global = new MKL_Complex16[num_trajectories * N];
	abs_rho_diag_omp = new double[num_omp_threads * N];
	rho_omp = new MKL_Complex16[num_omp_threads * N * N];
	etas = new double[num_trajectories];
	period = 0;

	for (int trajectoty_id = 0; trajectoty_id < num_trajectories; trajectoty_id++)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			phi_global[trajectoty_id * N + state_id].real = 0.0;
			phi_global[trajectoty_id * N + state_id].imag = 0.0;
		}

		etas[trajectoty_id] = 0.0;
	}

	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
		{
			abs_rho_diag_omp[thread_id * (N) + state_id_1] = 0.0;

			for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
			{
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].real = 0.0;
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].imag = 0.0;
			}
		}
	}
}

void refresh_dump_data(int N,
					   int num_omp_threads)
{
	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
		{
			abs_rho_diag_omp[thread_id * (N) + state_id_1] = 0.0;

			for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
			{
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].real = 0.0;
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].imag = 0.0;
			}
		}
	}
}

inline void recovery_phi_full(MKL_Complex16 * phi, double eta, double * g, void * A, int N, unsigned int k, VSLStreamStatePtr streamRand)  // MKL_complex16 * A
{
	int index = 0;
	MKL_Complex16 * res = new MKL_Complex16[N];
	double norm = sqrt(norm_vector2(phi, N));
	for(int i=0; i<N; i++)
	{
		phi[i].real /= (norm);
		phi[i].imag /= (norm);
	}
	norm=0.;

	double * gnorms = new double[k];
	double tmp = 0;
	double ran;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streamRand, 1, &ran, 0, 1);

	for(unsigned int i = 0; i < k; i++)
	{
		cblas_zgemv (CblasRowMajor, CblasNoTrans, N, N, &ONE, &((MKL_Complex16*)A)[IND2(i, 0)], N, phi, 1, &ZERO, res, 1);
		gnorms[i]=(norm_vector2(res, N));
		gnorms[i] *= g[i];
		tmp += gnorms[i];
	}

	ran *= tmp;

	while(ran - gnorms[index] > 0.)
	{
		ran -= gnorms[index];
		index++;
		if(index == k - 1)
			break;
	}

	while(gnorms[index] == 0)
	{
		if (index == 0)
		{
			index++;
		}
		else
		{
			index--;
		}
	}

	memset(res, 0, N*sizeof(MKL_Complex16));
	cblas_zgemv (CblasRowMajor, CblasNoTrans,N, N, &ONE, &(((MKL_Complex16*)A)[IND2(index,0)]), N, phi, 1, &ZERO, res, 1);

	for(int i=0; i<N; i++)
	{
		phi[i].real=res[i].real/sqrt(gnorms[index]/g[index]);
		phi[i].imag=res[i].imag/sqrt(gnorms[index]/g[index]);
	}

	delete (res);
	delete (gnorms);
}

inline void QJ_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, unsigned int N)
{
	cblas_zgemv (CblasRowMajor, CblasNoTrans, N, N, &ONE, matrix, N, phi, 1, &ZERO, res, 1);
}

void QJ_EXP_one_branch(MKL_Complex16 * phi, MKL_Complex16 * tmp_vec, double * eta, double * g, void * A, unsigned int k, split * branch,  VSLStreamStatePtr streamRand)
{
	if (branch->next == 0)
	{
		while (branch->counter != branch->steps)
		{
			QJ_step(phi, branch->matrix, tmp_vec, branch->N);
			if(if_norm(tmp_vec, eta, branch->N))
			{
				recovery_phi_full(tmp_vec, eta[0], g, A, branch->N, k, streamRand);
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streamRand, 1, eta, 0, 1);
				while (eta[0] == 0.0)
				{
					vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streamRand, 1, eta, 0, 1);
				}
			}

			memcpy(phi, tmp_vec, sizeof(MKL_Complex16)*branch->N);
			branch->counter++;

			t=t+branch->dt;
		}
	}
	else
	{
		while (branch->counter != branch->steps)
		{
			QJ_step(phi, branch->matrix, tmp_vec, branch->N);
			if(if_norm(tmp_vec, eta, branch->N))
			{
				QJ_EXP_one_branch(phi, tmp_vec, eta, g, A, k, branch->next, streamRand);
				t=t-branch->dt;
			}
			else
			{
				memcpy(phi, tmp_vec, sizeof(MKL_Complex16)*branch->N);
			}
			branch->counter++;
			t=t+branch->dt;
			if(max_steps == branch->steps)
			{
				unsigned int num_thread = omp_get_thread_num();
				unsigned int num_threads = omp_get_num_threads();
				int N = branch->N;

				double norm = norm_vector2(phi, N);
				MKL_Complex16 tmp;
				for(int k1 = 0; k1 < N; k1++)
				{
					for(int k2 = 0; k2 < N; k2++)
					{
						tmp = Complex_mul(Complex_scalar_mul(&phi[k1], &phi[k2], 1), 1/norm);
						matrix_omp[(index_kT[num_thread] + num_thread * gkT * 32) * N * N + k1 * N + k2].real += tmp.real;
						matrix_omp[(index_kT[num_thread] + num_thread * gkT * 32) * N * N + k1 * N + k2].imag += tmp.imag;


						matrix_omp[(index_kT[num_thread] + num_thread * gkT * 32 + num_threads * gkT * 32) * N * N + k1 * N + k2].real += tmp.real * tmp.real;
						matrix_omp[(index_kT[num_thread] + num_thread * gkT * 32 + num_threads * gkT * 32) * N * N + k1 * N + k2].imag += tmp.imag * tmp.imag;
					}
				}
				index_kT[num_thread]++;
			}


		}
	}

	branch->counter = 0;
}

void qj_propagate_one_period(MKL_Complex16 * phi, double * eta, split * head,  VSLStreamStatePtr streamRand)
{
	MKL_Complex16 * tmp_vec = new MKL_Complex16[head->N];
	for (unsigned int i = 0; i < head->counter; i++)
	{
		QJ_EXP_one_branch(phi, tmp_vec, eta, head->g, head->matrix, head->steps, &(head->next)[i], streamRand);
	}
	delete(tmp_vec);
}

void set_init_conditions(MKL_Complex16 * phi, int N, int state_id)
{
	if (state_id >= 0 && state_id < N)
	{
		phi[state_id].real = 1.0;
	}
	else
	{
		for (int i = 0; i < N; i++)
		{
			phi[i].real = sqrt(1.0/double(N));
		}
	}
}

void calc_and_dump_characteristics(int N, int Nc, char* write_type)
{
	char characteristics_file_name[] = "characteristics.txt";

	double start_time = omp_get_wtime();

	MKL_Complex16 * tmp_mult = new MKL_Complex16[N * N];

	for(int state_1 = 0; state_1 < N; state_1++)
	{
		for(int state_2 = 0; state_2 < N; state_2++)
		{
			tmp_mult[state_1 * N + state_2].real = 0.0;
			tmp_mult[state_1 * N + state_2].imag = 0.0;
		}
	}

	MKL_Complex16 alpha;
	alpha.real = 1.0;
	alpha.imag = 0.0;

	MKL_Complex16 beta;
	beta.real = 0.0;
	beta.imag = 0.0;

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &alpha,
		ev_t_hamiltonian, N, rho_curr, N, &beta, tmp_mult, N);

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &alpha,
		tmp_mult, N, ev_hamiltonian, N, &beta, rho_curr_and, N);

	delete[] tmp_mult;

	// ================== Entropy ========================
	MatrixXcd matrix_rho_eigen(N, N);
	MatrixXcd matrix_rho_and_eigen(N, N);
	for(int state_1 = 0; state_1 < N; state_1++)
	{
		for(int state_2 = 0; state_2 < N; state_2++)
		{
			std::complex<double> val;
			val.real(rho_curr[state_1 * N + state_2].real);
			val.imag(rho_curr[state_1 * N + state_2].imag);
			matrix_rho_eigen(state_1, state_2) = val;

			std::complex<double> val_and;
			val_and.real(rho_curr_and[state_1 * N + state_2].real);
			val_and.imag(rho_curr_and[state_1 * N + state_2].imag);
			matrix_rho_and_eigen(state_1, state_2) = val_and;
		}
	}

	MatrixXcd matrix_log_rho_eigen = matrix_rho_eigen.log();
	MatrixXcd matrix_log_rho_and_eigen = matrix_rho_and_eigen.log();

	MKL_Complex16 * matrix_log_rho_mkl = new MKL_Complex16[N*N];
	MKL_Complex16 * mult_result = new MKL_Complex16[N*N];
	MKL_Complex16 * matrix_log_rho_and_mkl = new MKL_Complex16[N*N];
	MKL_Complex16 * mult_result_and = new MKL_Complex16[N*N];
	for(int state_1 = 0; state_1 < N; state_1++)
	{
		for(int state_2 = 0; state_2 < N; state_2++)
		{
			std::complex<double> val = matrix_log_rho_eigen(state_1, state_2);
			matrix_log_rho_mkl[state_1 * N + state_2].real = val.real();
			matrix_log_rho_mkl[state_1 * N + state_2].imag = val.imag();

			mult_result[state_1 * N + state_2].real = 0.0;
			mult_result[state_1 * N + state_2].imag = 0.0;

			std::complex<double> val_and = matrix_log_rho_and_eigen(state_1, state_2);
			matrix_log_rho_and_mkl[state_1 * N + state_2].real = val_and.real();
			matrix_log_rho_and_mkl[state_1 * N + state_2].imag = val_and.imag();

			mult_result_and[state_1 * N + state_2].real = 0.0;
			mult_result_and[state_1 * N + state_2].imag = 0.0;
		}
	}

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, rho_curr, N, matrix_log_rho_mkl, N, &ZERO, mult_result, N);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, rho_curr_and, N, matrix_log_rho_and_mkl, N, &ZERO, mult_result_and, N);

	entropy = 0.0;
	entropy_and = 0.0;
	for(int state = 0; state < N; state++)
	{
		entropy -=  mult_result[state * N + state].real;
		entropy_and -= mult_result_and[state * N + state].real;
	}

	delete[] matrix_log_rho_mkl;
	delete[] mult_result;
	delete[] matrix_log_rho_and_mkl;
	delete[] mult_result_and;
	// ====================================================

	// ============= Imbalance ==================
	double * n_part = new double[Nc];
	double * n_part_and = new double[Nc];

	for (int i = 0; i < Nc; i++)
	{
		n_part[i] = 0.0;
		n_part_and[i] = 0.0;
	}

	for (int state_id = 0; state_id < N; state_id++)
	{
		int state_number = state_in_bin[state_id];

		for (int i = 0; i < Nc; i++)
		{
			int state_bin_code = state_number & (1 << i) ? 1 : 0;
			n_part[i] += rho_curr[state_id * N + state_id].real * double(state_bin_code);
			n_part_and[i] += rho_curr_and[state_id * N + state_id].real * double(state_bin_code);

		}
	}

	double sum_odd = 0.0;
	double sum_even = 0.0;
	double sum_all = 0.0;

	double sum_odd_and = 0.0;
	double sum_even_and = 0.0;
	double sum_all_and = 0.0;

	for (int i = 0; i < Nc; i++)
	{
		if (i%2 == 0)
		{
			sum_odd += n_part[i];
			sum_odd_and += n_part_and[i];
		}
		else
		{
			sum_even += n_part[i];
			sum_even_and += n_part_and[i];
		}
	}
	sum_all = sum_even + sum_odd;
	sum_all_and = sum_even_and + sum_odd_and;

	imbalance = (sum_odd - sum_even) / sum_all;
	imbalance_and = (sum_odd_and - sum_even_and) / sum_all_and;

	delete[] n_part;
	delete[] n_part_and;
	// ==============================================

	// =================== IPR ======================
	ipr = 0.0;
	ipr_and = 0.0;
	for (int state_id = 0; state_id < N; state_id++)
	{
		// diagonal elements contains only real part
		ipr += (rho_curr[state_id * N + state_id].real * rho_curr[state_id * N + state_id].real + rho_curr[state_id * N + state_id].imag * rho_curr[state_id * N + state_id].imag);
		ipr_and += (rho_curr_and[state_id * N + state_id].real * rho_curr_and[state_id * N + state_id].real + rho_curr_and[state_id * N + state_id].imag * rho_curr_and[state_id * N + state_id].imag);
	}
	// ==============================================

	FILE * characteristics_file = fopen(characteristics_file_name, write_type);
	fprintf(characteristics_file, "%0.18le %0.18le %0.18le %0.18le %0.18le %0.18le\n", entropy, imbalance, ipr, entropy_and, imbalance_and, ipr_and);
	fclose(characteristics_file);

	double time = omp_get_wtime() - start_time;
	printf("routines: %0.3le\n", time);
}

void omp_qj_single_trajectory_init_propagation(split * head,
											   VSLStreamStatePtr rnd_stream,
											   MKL_Complex16 * phi,
											   double * abs_rho_diag,
											   MKL_Complex16 * rho,
											   int num_periods_in_trans_proc,
											   int num_trajectories,
											   int trajectory_id,
											   int init_state_id,
											   int thread_id)
{
	int N = head->N;

	set_init_conditions(phi, N, init_state_id);

	double eta = 0.0;
	while (eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd_stream, 1, &eta, 0, 1);
	}

	if (num_periods_in_trans_proc > 0)
	{
		for(int period_id = 0; period_id < num_periods_in_trans_proc; period_id++)
		{
			qj_propagate_one_period(phi, &eta, head, rnd_stream);
		}
	}
	else
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream);
	}

	etas[trajectory_id] = eta;

	double norm = norm_vector2(phi, N);
	for(int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		abs_rho_diag[state_id_1] += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_1], 1), 1.0 / norm).real / num_trajectories;

		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho[state_id_1 * (N) + state_id_2].real += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).real / num_trajectories;
			rho[state_id_1 * (N) + state_id_2].imag += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).imag / num_trajectories;
		}
	}
}

void dump_propagation_data(int N,
						   int Nc,
						   int num_omp_threads,
						   char* write_type,
						   int calc_characteristics,
						   int dump_rho)
{
	char abs_rho_diag_file_name[] = "abs_rho_diag.txt";
	char rho_file_name[512];
	sprintf(rho_file_name, "rho_period_%d.txt", period);
	char periods_file_name[] = "periods.txt";

	FILE * abs_rho_diag_file = fopen(abs_rho_diag_file_name, write_type);
	double norm = 0.0;
	for (int state_id = 0; state_id < N; state_id++)
	{
		double abs_rho_diag_val = 0.0;
		for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
		{
			abs_rho_diag_val += abs_rho_diag_omp[thread_id * (N) + state_id];
		}
		fprintf(abs_rho_diag_file, "%0.18le\n", abs_rho_diag_val);
		norm += abs_rho_diag_val;
	}
	printf("diff norm: %0.18le\n", 1.0 - norm);
	fclose(abs_rho_diag_file);

	FILE * rho_file;
	if (dump_rho)
	{
		rho_file = fopen(rho_file_name, "w");
	}
	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			double real_part = 0.0;
			double imag_part = 0.0;
			for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
			{
				real_part += rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].real;
				imag_part += rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].imag;
			}

			rho_curr[state_id_1 * (N) + state_id_2].real = real_part;
			rho_curr[state_id_1 * (N) + state_id_2].imag = imag_part;

			if (dump_rho)
			{
				fprintf(rho_file, "%0.18le %0.18le\n", real_part, imag_part);
			}
		}
	}
	if (dump_rho)
	{
		fclose(rho_file);
	}

	FILE * periods_file = fopen(periods_file_name, write_type);
	fprintf(periods_file, "%d\n", period);
	fclose(periods_file);

	if (calc_characteristics)
	{
		calc_and_dump_characteristics(N, Nc, write_type);
	}

}






void omp_qj_single_trajectory_propagation_to_dump(split * head,
												  VSLStreamStatePtr rnd_stream,
												  MKL_Complex16 * phi,
												  double * abs_rho_diag,
												  MKL_Complex16 * rho,
												  int num_trajectories,
												  int trajectory_id,
												  int begin_propagation_loop,
												  int end_propagation_loop,
												  int thread_id)
{
	int N = head->N;

	double eta = etas[trajectory_id];

	for(int period_id = begin_propagation_loop; period_id < end_propagation_loop; period_id++)
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream);
	}

	etas[trajectory_id] = eta;

	double norm = norm_vector2(phi, N);
	for(int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		abs_rho_diag[state_id_1] += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_1], 1), 1.0 / norm).real / num_trajectories;

		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho[state_id_1 * (N) + state_id_2].real += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).real / num_trajectories;
			rho[state_id_1 * (N) + state_id_2].imag += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).imag / num_trajectories;
		}
	}
}

void init_aux_data(char aux_file_name[], int N)
{
	state_in_bin = new int[N];
	rho_curr = new MKL_Complex16[N*N];
	rho_curr_and = new MKL_Complex16[N*N];
	hamiltonian = new MKL_Complex16[N*N];
	ev_hamiltonian = new MKL_Complex16[N*N];
	ev_t_hamiltonian = new MKL_Complex16[N*N];

	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho_curr[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr[state_id_1 * (N) + state_id_2].imag = 0.0;

			rho_curr_and[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr_and[state_id_1 * (N) + state_id_2].imag = 0.0;

			hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;
		}
	}


	FILE * aux_file = fopen(aux_file_name, "rb");
	fread(state_in_bin, sizeof(int), N, aux_file);
	fread(hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(ev_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(ev_t_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fclose(aux_file);
}

void refresh_aux_data(int N)
{
	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho_curr[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr[state_id_1 * (N) + state_id_2].imag = 0.0;

			rho_curr_and[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr_and[state_id_1 * (N) + state_id_2].imag = 0.0;
		}
	}
}

void delete_aux_data()
{
	delete[] state_in_bin;
	delete[] rho_curr;
	delete[] rho_curr_and;
	delete[] hamiltonian;
	delete[] ev_hamiltonian;
	delete[] ev_t_hamiltonian;
}

void omp_qj_many_trajectories(char input_file_name[],
							  char aux_file_name[],
							  int Nc,
							  int num_periods,
							  int num_dumps,
							  int dump_type,
							  int num_periods_in_trans_proc,
							  int num_omp_threads,
							  int num_trajectories,
							  int rnd_max,
							  int rnd_cur,
							  int init_state_id,
							  int calc_characteristics,
							  int dump_rho)
{
	printf("num_omp_threads: %d\n\n", num_omp_threads);
	omp_set_num_threads(num_omp_threads);

	FILE * input_file;
	input_file = fopen(input_file_name, "rb");
	split * head = create_struct_bin(input_file);
	fclose(input_file);

	int N = head->N;

	VSLStreamStatePtr * rnd_streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream (&rnd_streams[0], VSL_BRNG_MCG59, 777);
	for(int i = 1; i < num_trajectories; i++)
	{
		vslCopyStream(&rnd_streams[i], rnd_streams[0]);
	}
	for(int i = 0; i < num_trajectories; i++)
	{
		vslLeapfrogStream(rnd_streams[i], rnd_cur + i, rnd_max);
	}

	split * heads = new split[num_omp_threads];
	for(int i = 0; i < num_omp_threads; i++)
	{
		cmp_struct_not_member(head, &heads[i]);
	}

	init_data(N, num_trajectories, num_omp_threads);
	init_aux_data(aux_file_name, N);

	double step = double(num_periods - num_periods_in_trans_proc) / double(num_dumps - 1);
	double start = double(num_periods_in_trans_proc);
	if (dump_type == 1)
	{
		start = 0.0;
		if (num_periods_in_trans_proc > 0)
		{
			start = log10(double(num_periods_in_trans_proc));
		}

		step = (log10(double(num_periods)) - start) / double(num_dumps - 1);
	}

#pragma omp parallel for
	for(int i = 0; i < num_trajectories; i++)
	{
		int thread_id = omp_get_thread_num();

		omp_qj_single_trajectory_init_propagation(
			&heads[thread_id],
			rnd_streams[i],
			&phi_global[i * N],
			&abs_rho_diag_omp[thread_id * N],
			&rho_omp[thread_id * N * N],
			num_periods_in_trans_proc,
			num_trajectories,
			i,
			init_state_id,
			thread_id
			);
	}

	if (num_periods_in_trans_proc > 0)
	{
		period = num_periods_in_trans_proc;
	}
	else
	{
		period = 1;
	}
	printf("Period = %d\n", period);

	dump_propagation_data(N, Nc, num_omp_threads, "w", calc_characteristics, dump_rho);
	refresh_dump_data(N, num_omp_threads);
	refresh_aux_data(N);

	double current_dump_limit = start;
	int begin_propagation_loop = period;
	int end_propagation_loop = 0;

	while (begin_propagation_loop != num_periods)
	{
		current_dump_limit += step;
		end_propagation_loop = num_periods;

		if (dump_type == 0)
		{
			if (current_dump_limit < double(end_propagation_loop))
			{
				end_propagation_loop = int(current_dump_limit);
			}
		}
		else if (dump_type == 1)
		{
			int exp_limit = int(pow(10.0, current_dump_limit + 1.0e-10));
			if (exp_limit < end_propagation_loop)
			{
				if (exp_limit <= begin_propagation_loop)
				{
					end_propagation_loop = begin_propagation_loop + 1;
				}
				else
				{
					end_propagation_loop = exp_limit;
				}
			}
		}

#pragma omp parallel for
		for(int i = 0; i < num_trajectories; i++)
		{
			int thread_id = omp_get_thread_num();

			omp_qj_single_trajectory_propagation_to_dump(
				&heads[thread_id],
				rnd_streams[i],
				&phi_global[i * N],
				&abs_rho_diag_omp[thread_id * N],
				&rho_omp[thread_id * N * N],
				num_trajectories,
				i,
				begin_propagation_loop,
				end_propagation_loop,
				thread_id
				);
		}

		printf("Period = %d\n", end_propagation_loop);
		begin_propagation_loop = end_propagation_loop;
		period = end_propagation_loop;

		dump_propagation_data(N, Nc, num_omp_threads, "a", calc_characteristics, dump_rho);
		refresh_dump_data(N, num_omp_threads);
		refresh_aux_data(N);
	}

	delete_aux_data();

	for(int i = 0; i < num_omp_threads; i++)
	{
		delete_split_struct_not_member (&heads[i]);
	}
	delete (heads);
	delete_split_struct (head);
	delete (rnd_streams);
	delete[] phi_global;
	delete[] abs_rho_diag_omp;
	delete[] rho_omp;
	delete[] etas;
}
