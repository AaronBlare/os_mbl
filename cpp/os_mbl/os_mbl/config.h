#pragma once
#include "config.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <complex>
#include <mkl.h>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <math.h>

using namespace std;

struct ConfigParam
{
	int   task;						// task type: 0 - zev; 1 - int; 2 - tr; 3 - qj
	int		Nc;						// num of cells
	int		dt;						// dissipator type
	double	dp;						// Diehl dissipator phase
	int		et;						// energy type
	int		bc;						// border conditions
	double	W;						// disorder
	double	U;						// interaction
	double	J;						// hopping
	double	g;						// dissipation
	int		seed;					// seed

	int		max_num_seeds;			// max num of seeds
	int		dump_vecs;				// dump vectors?
	int		dump_mtxs;				// dump matrixes?
	int		dump_sops;				// dump superoperators?
	string	init_file;					// energies init file (if nessesary)
	string  path;					// dump path (if nessesary)

	ConfigParam(
		int _task				= 0,
		int	_Nc					= 10,
		int	_dt					= 1,
		double _dp				= 0.0,
		int	_et					= 0,
		int _bc					= 0,
		double	_W				= 8.0,
		double	_U				= 1.0,
		double	_J				= 1.0,
		double	_g				= 0.1,	
		int		_seed			= 1,

		int _max_num_seeds		= 10000,
		int	_dump_vecs			= 0,
		int _dump_mtxs			= 0,
		int _dump_sops			= 0,
		string _init_file		= "",
		string _path			= ""
		)
	{
		task						= _task;
		Nc							= _Nc;
		dt							= _dt;
		dp							= _dp;
		et							= _et;
		bc							= _bc;
		W							= _W;
		U							= _U;
		J							= _J;
		g							= _g;
		seed						= _seed;

		max_num_seeds				= _max_num_seeds;
		dump_vecs					= _dump_vecs;
		dump_mtxs					= _dump_mtxs;
		dump_sops					= _dump_sops;
		init_file					= _init_file;
		path						= _path;
	}
};


void set_param(ConfigParam &param, string str, string value);

void init_config_param(ConfigParam &param, char * file_name);

void output_setting(ConfigParam &param);
