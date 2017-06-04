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

	if (param.task == 0)
	{
		run_zero_eigen_vector(cd, param);
	}
	else if (param.task == 1)
	{
		run_f_ode(cd, param);
	}
	else if (param.task == 2)
	{
		run_trans_rates(cd, param);
	}
	else
	{
		stringstream msg;
		msg << "wrong task value: " << param.task << endl;
		Error(msg.str());
	}

	time = omp_get_wtime() - time;
	cout << "total time: " << time << endl << endl;

	return 0;
}