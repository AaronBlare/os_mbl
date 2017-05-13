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
#include <boost/numeric/odeint.hpp>

using namespace std;

typedef vector< complex < double > > state_type;

using namespace boost::numeric::odeint;

struct ConfigParam
{
	int   task;						// task type: 0 - zev; 1 - int; 2 - tr; 3 - qj
	int		Nc;						// num of cells
	double	W;						// disorder
	double	U;						// interaction
	double	J;						// hopping
	double	g;						// dissipation
	int		dt;						// dissipator type
	double	alpha;					// Diehl dissipator phase
	int		et;						// energy type
	int		seed;					// random seed
	int		max_num_seeds;			// max num of seeds
	int		dump_vecs;				// dump vectors?
	int		dump_mtxs;				// dump matrixes?
	int		dump_sops;				// dump superoperators?
	string	init;					// energies init file (if nessesary)
	string  path;					// dump path (if nessesary)

	int init_state_type;			// integration: 0 - localized in init_state_id; 1 - uniform
	int init_state_id;				// integration: init state localization index (if init_state_type == 0)
	int int_dump_type;				// integration: 0 - linear dumps; 1 - log dumps
	int num_dumps;					// integration: number of dumps
	double begin_dump_time;			// integration: first dumped time
	double end_dump_time;			// integration: last dumped time

	ConfigParam(
		int _task						= 0,
		int		_Nc						= 10,
		double	_W						= 8.0,
		double	_U						= 1.0,
		double	_J						= 1.0,
		double	_g						= 0.1,
		int		_dt						= 1,
		double _alpha					= 0.0,
		int		_et						= 0,
		int		_seed					= 1,
		int _max_num_seeds				= 10000,
		int	_dump_vecs					= 0,
		int _dump_mtxs					= 0,
		int _dump_sops					= 0,
		string _init					= "",
		string _path					= "",
		int _init_state_type			= 0,			
		int _init_state_id				= 49,				
		int _int_dump_type				= 1,
		int _num_dumps					= 50,					
		double _begin_dump_time			= 0.1,			
		double _end_dump_time			= 10000.0
		)
	{
		task						= _task;
		Nc							= _Nc;
		W							= _W;
		U							= _U;
		J							= _J;
		g							= _g;
		dt							= _dt;
		alpha						= _alpha;
		et							= _et;
		seed						= _seed;
		max_num_seeds				= _max_num_seeds;
		dump_vecs					= _dump_vecs;
		dump_mtxs					= _dump_mtxs;
		dump_sops					= _dump_sops;
		init						= _init;
		path						= _path;

		init_state_type = _init_state_type;
		init_state_id = _init_state_id;
		int_dump_type = _int_dump_type;
		num_dumps = _num_dumps;
		begin_dump_time = _begin_dump_time;
		end_dump_time = _end_dump_time;
	}
};


void set_param(ConfigParam &param, string str, string value);

void init_config_param(ConfigParam &param, char * file_name);

void output_setting(ConfigParam &param);
