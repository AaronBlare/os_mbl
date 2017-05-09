#include "mkl.h"
#include "stdio.h"
#include "math.h"
#include "string.h"


struct split
{
	bool type;
	split * next;
	split * prev;
	double dt;
	unsigned int steps;
	unsigned int counter;
	unsigned int N;
	MKL_Complex16 * matrix;
	double * g;
};

split * create_struct_bin (FILE * file);
split * create_struct (FILE * file);

void delete_split_struct (split * head);
void qj_propagate_one_period(MKL_Complex16 * phi, double eta, double * g, split * head,  VSLStreamStatePtr streamRand);
double norm_vector2(MKL_Complex16 * phi, int N);

void cmp_struct_not_member(split * head1, split * head2);
void delete_split_struct_not_member (split * head);

void set_init_conditions(MKL_Complex16 * phi, int N, int state_id);

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
							  int dump_rho);

