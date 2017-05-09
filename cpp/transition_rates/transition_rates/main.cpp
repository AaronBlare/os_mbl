#include <omp.h>
#include "config.h"
#include "data.h"


int main(int argc, char ** argv)
{
	cout << "current path: " << argv[0] << endl << endl;
	if (argc > 1)
	{
		omp_set_num_threads(atoi(argv[1]));
	}
	ConfigParam param;
	init_config_param(param, "config.txt");

	ConfigData cd;

	double time = omp_get_wtime();

	run_experiment(cd, param);

	time = omp_get_wtime() - time;
	cout << "total time: " << time << endl << endl;

	return 0;
}