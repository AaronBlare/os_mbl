#include "utils.h"

void error(const string& err, const char* func, const char* file, int line)
{
	cout << "Error: " << err << " in " << func << " (" << file << ", line: " << line << ")" << endl;
	exit(EXIT_FAILURE);
}

int n_choose_k(int n, int k)
{
	int res = 1;

	if ( k > n - k )
	{
		k = n - k;
	}

	for (int i = 0; i < k; ++i)
	{
		res *= (n - i);
		res /= (i + 1);
	}

	return res;
}

int bit_count(int value) 
{
	int count = 0;

	while (value > 0) 
	{
		if ((value & 1) == 1)
		{
			count++;
		}
		value >>= 1;
	}

	return count;
}

int bit_at(int value, int position)
{
	if ((value >> position) & 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void print_int_array(int * data, int N)
{
	cout << endl;
	for (int i = 0; i < N; i ++)
	{
		cout << i << ": " << data[i] << endl;
	}
}

string file_name_suffix(ConfigParam &cp, int precision)
{
	if (cp.task == 0)
	{
		return file_name_suffix_zev(cp, precision);
	}
	else if (cp.task == 1)
	{
		return file_name_suffix_f_ode(cp, precision);
	}
	else if (cp.task == 2)
	{
		return file_name_suffix_f_int(cp, precision);
	}
	else if (cp.task == 3)
	{
		return file_name_suffix_tr(cp, precision);
	}
	else
	{
		stringstream msg;
		msg << "wrong task value: " << cp.task << endl;
		Error(msg.str());
	}
}

string file_name_suffix_zev(ConfigParam &cp, int precision)
{
	stringstream fns;
	fns << "_zev";
	fns << "_Nc(" << cp.Nc << ")";
	fns << "_dt(" << cp.dt << ")";
	fns << "_dp(" << setprecision(precision) << fixed << cp.dp << ")";
	fns << "_et(" << cp.et << ")";
	fns << "_bc(" << cp.bc << ")";
	fns << "_W(" << setprecision(precision) << fixed << cp.W << ")";
	fns << "_U(" << setprecision(precision) << fixed << cp.U << ")";
	fns << "_J(" << setprecision(precision) << fixed << cp.J << ")";
	fns << "_g(" << setprecision(precision) << fixed << cp.g << ")";
	fns << "_seed(" << cp.seed << ")";
	fns << ".txt";

	return fns.str();
}

string file_name_suffix_f_ode(ConfigParam &cp, int precision)
{
	stringstream fns;
	fns << "_f_ode";
	fns << "_Nc(" << cp.Nc << ")";
	fns << "_dt(" << cp.dt << ")";
	fns << "_dp(" << setprecision(precision) << fixed << cp.dp << ")";
	fns << "_et(" << cp.et << ")";
	fns << "_bc(" << cp.bc << ")";
	fns << "_W(" << setprecision(precision) << fixed << cp.W << ")";
	fns << "_U(" << setprecision(precision) << fixed << cp.U << ")";
	fns << "_J(" << setprecision(precision) << fixed << cp.J << ")";
	fns << "_g(" << setprecision(precision) << fixed << cp.g << ")";
	fns << "_seed(" << cp.seed << ")";
	fns << ".txt";

	return fns.str();
}

string file_name_suffix_f_int(ConfigParam &cp, int precision)
{
	stringstream fns;
	fns << "_f_int";
	fns << "_Nc(" << cp.Nc << ")";
	fns << "_dt(" << cp.dt << ")";
	fns << "_dp(" << setprecision(precision) << fixed << cp.dp << ")";
	fns << "_et(" << cp.et << ")";
	fns << "_bc(" << cp.bc << ")";
	fns << "_W(" << setprecision(precision) << fixed << cp.W << ")";
	fns << "_U(" << setprecision(precision) << fixed << cp.U << ")";
	fns << "_J(" << setprecision(precision) << fixed << cp.J << ")";
	fns << "_g(" << setprecision(precision) << fixed << cp.g << ")";
	fns << "_seed(" << cp.seed << ")";
	fns << ".txt";

	return fns.str();
}

string file_name_suffix_tr(ConfigParam &cp, int precision)
{
	stringstream fns;
	fns << "_tr";
	fns << "_Nc(" << cp.Nc << ")";
	fns << "_dt(" << cp.dt << ")";
	fns << "_dp(" << setprecision(precision) << fixed << cp.dp << ")";
	fns << "_et(" << cp.et << ")";
	fns << "_bc(" << cp.bc << ")";
	fns << "_W(" << setprecision(precision) << fixed << cp.W << ")";
	fns << "_U(" << setprecision(precision) << fixed << cp.U << ")";
	fns << "_J(" << setprecision(precision) << fixed << cp.J << ")";
	fns << "_g(" << setprecision(precision) << fixed << cp.g << ")";
	fns << "_seed(" << cp.seed << ")";
	fns << ".txt";

	return fns.str();
}

void write_double_data(string file_name, double * data, int size, int precision, bool append)
{
	if (append)
	{
		ofstream ofs(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else 
	{
		ofstream ofs(file_name);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

void write_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append)
{
	if (append)
	{
		ofstream ofs(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i].real << " " << data[i].imag << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else
	{
		ofstream ofs(file_name);
		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i].real << " " << data[i].imag << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

void write_sparse_complex_mtx(string file_name, crsMatrix *A, int precision, bool append)
{
	if (append)
	{
		ofstream ofs(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int i = 0; i < A->N; i++)
			{
				for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
				{
					ofs << i + 1 << " " << A->Col[k] + 1 << " " << A->Value[k].re << " " << A->Value[k].im << endl;
				}
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else
	{
		ofstream ofs(file_name);
		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int i = 0; i < A->N; i++)
			{
				for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
				{
					ofs << i + 1 << " " << A->Col[k] + 1 << " " << A->Value[k].re << " " << A->Value[k].im << endl;
				}
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

vector<int> convert_int_to_vector_of_bits(int x, int size) 
{
	vector<int> res;

	int id = 0;

	while(id < size) 
	{
		if (x&1)
		{
			res.push_back(1);
		}
		else
		{
			res.push_back(0);
		}

		x >>= 1;
		id ++;
	}

	reverse(res.begin(), res.end());

	return res;
}

vector<int> sort_doubles_with_order(vector<double> &v)
{
	vector<int> order(v.size());

	iota(order.begin(), order.end(), 0);

	sort(order.begin(), order.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return order;
}