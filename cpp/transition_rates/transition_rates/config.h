#pragma once
#include "config.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <mkl.h>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <math.h>

using namespace std;

struct ConfigParam
{
	int		Nc;				// num of cells
	double	W;				// disorder
	double	U;				// interaction
	double	J;				// hopping
	double	g;				// dissipation
	int		dt;				// dissipator type
	double	alpha;			// Diehl dissipator phase
	int		et;				// energy type
	string	init;			// random init file
	int		seed;			// random seed
	int		max_num_seeds;	// max num of seeds
	int		dump;			// dump type
	string  path;

	ConfigParam(
		int		_Nc			= 10,
		double	_W			= 8.0,
		double	_U			= 1.0,
		double	_J			= 1.0,
		double	_g			= 0.1,
		int		_dt			= 1,
		double _alpha		= 0.0,
		int		_et			= 0,
		string _init		= "",
		int		_seed		= 1,
		int _max_num_seeds	= 10000,
		int		_dump		= 0,
		string _path		= ""
		)
	{
		Nc			= _Nc;
		W			= _W;
		U			= _U;
		J			= _J;
		g			= _g;
		dt			= _dt;
		alpha		= _alpha;
		et			= _et;
		init		= _init;
		seed		= _seed;
		max_num_seeds	= _max_num_seeds;
		dump		= _dump;
		path = _path;
	}
};


void set_param(ConfigParam &param, string str, string value);

void init_config_param(ConfigParam &param, char * file_name);

void output_setting(ConfigParam &param);
