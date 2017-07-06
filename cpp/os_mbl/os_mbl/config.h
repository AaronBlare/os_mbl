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

	double		int_h;   // step
	int			int_ist; // init state type
	int			int_isi; // init state index
	int			int_dt;	 // dump type
	double		int_db;  // begin dump
	double		int_de;	 // end dump
	int			int_dn;	 // num dumps

	int		max_num_seeds;			// max num of seeds
	int		dump_vecs;				// dump vectors?
	int		dump_mtxs;				// dump matrixes?
	int		dump_sops;				// dump superoperators?
	string	init_file;				// energies init file (if nessesary)
	string  path;					// dump path (if nessesary)

	ConfigParam(
		int _task = 0,
		int	_Nc = 10,
		int	_dt = 1,
		double _dp = 0.0,
		int	_et = 0,
		int _bc = 0,
		double	_W = 8.0,
		double	_U = 1.0,
		double	_J = 1.0,
		double	_g = 0.1,
		int		_seed = 1,

		double	_int_h = 0.01,
		int		_int_ist = 0,
		int		_int_isi = 0,
		int		_int_dt = 1,
		double	_int_db = 0,
		double	_int_de = 1,
		int		_int_dn = 1,

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
		
		int_h = _int_h;
		int_ist = _int_ist;
		int_isi = _int_isi;
		int_dt = _int_dt;
		int_db = _int_db;
		int_de = _int_de;
		int_dn = _int_dn;

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
