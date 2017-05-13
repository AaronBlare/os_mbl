#pragma once
#include "config.h"
#include "utils.h"

class ode_lindbladian 
{
	int dim;

	MKL_Complex16 * lindbladian;
	MKL_Complex16 * curr;
	MKL_Complex16 * next;
	
	MKL_Complex16 ONE;
	MKL_Complex16 ZERO;

public:
	ode_lindbladian(int _dim, MKL_Complex16 * _lindbladian, MKL_Complex16 * _curr, MKL_Complex16 * _next);

	void operator() (const state_type &x, state_type &dxdt, const double t);

};

struct observer
{	
	double t_curr;
	double t_next;
	vector<double> & times;

	observer(std::vector< double > &t)
		: times(t) { }

	void operator()(const state_type &x, double t);
};