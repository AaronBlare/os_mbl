#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include "header.h"
#include "omp.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		printf("Not enough arguments\n");
		return -2;
	}

	int Nc;
	int num_periods;
	int init_state_id;
	int num_dumps;
	int dump_type;
	int num_periods_in_trans_proc;
	int num_omp_threads;
	int num_trajectories;
	int rnd_max;
	int rnd_cur;
	int calc_characteristics;
	int dump_rho;
	
	FILE * config_file = fopen("config.txt", "r");
	fscanf(config_file, "Nc = %d\n", &Nc);
	fscanf(config_file, "num_periods = %d\n", &num_periods);
	fscanf(config_file, "init_state_id = %d\n", &init_state_id);
	fscanf(config_file, "num_dumps = %d\n", &num_dumps);
	fscanf(config_file, "dump_type = %d\n", &dump_type);
	fscanf(config_file, "num_periods_in_trans_proc = %d\n", &num_periods_in_trans_proc);
	fscanf(config_file, "num_omp_threads = %d\n", &num_omp_threads);
	fscanf(config_file, "num_trajectories = %d\n", &num_trajectories);
	fscanf(config_file, "rnd_max = %d\n", &rnd_max);
	fscanf(config_file, "rnd_cur = %d\n", &rnd_cur);
	fscanf(config_file, "calc_characteristics = %d\n", &calc_characteristics);
	fscanf(config_file, "dump_rho = %d\n", &dump_rho);
	fclose (config_file);

	double start_time = omp_get_wtime();
	omp_qj_many_trajectories(argv[1], argv[2], Nc, num_periods, num_dumps, dump_type, num_periods_in_trans_proc, num_omp_threads, num_trajectories, rnd_max, rnd_cur, init_state_id, calc_characteristics, dump_rho);
	double time = omp_get_wtime() - start_time;

	printf("elapsed_time: %0.3le\n", time);

	return 0;
}