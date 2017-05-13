#include "integration.h"

ode_lindbladian::ode_lindbladian(int _dim, MKL_Complex16 * _lindbladian, MKL_Complex16 * _curr, MKL_Complex16 * _next)
{
	dim = _dim;

	lindbladian = _lindbladian;
	curr = _curr;
	next = _next;

	ONE = { 1.0, 0.0 };
	ZERO = { 0.0, 0.0 };
}

void ode_lindbladian::operator() (const state_type &x, state_type &dxdt, const double t)
{
	for (int state_id = 0; state_id < dim; state_id++)
	{
		curr[state_id].real = x[state_id].real();
		curr[state_id].imag = x[state_id].imag();

		next[state_id].real = 0.0;
		next[state_id].imag = 0.0;
	}

	cblas_zgemv(
		CblasRowMajor,
		CblasNoTrans,
		dim,
		dim,
		&ONE,
		lindbladian,
		dim,
		curr,
		1,
		&ZERO,
		next,
		1
	);

	for (int state_id = 0; state_id < dim; state_id++)
	{
		dxdt[state_id].real(next[state_id].real);
		dxdt[state_id].imag(next[state_id].imag);
	}
}

void observer::operator()(const state_type &x, double t)
{
	times.push_back(t);
}
