#include "data.h"

void run_experiment(ConfigData &cd, ConfigParam &cp)
{
	cd.Nc = cp.Nc;

	if (cd.Nc % 2 != 0)
	{
		stringstream msg;
		msg << "Nc must divide by 2 without any remainder" << endl;
		Error(msg.str());
	}

	cd.Np = cd.Nc / 2;
	cd.Ns = n_choose_k(cd.Nc, cd.Np);

	cout << "num states = " << cd.Ns << endl;

	int num_states = init_aux_data(cd);

	if (cd.Ns != num_states)
	{
		stringstream msg;
		msg << "wrong num states calculation" << endl;
		Error(msg.str());
	}

	init_hamiltonian_data(cd, cp);
	init_hamiltonian_eig_data(cd, cp);
	free_hamiltonian_data(cd);
	init_trans_rates_data(cd, cp);
	init_diag_rho_data(cd, cp);
	free_hamiltonian_eig_data(cd);
	free_trans_rates_data(cd);

	calculate_characteristics(cd, cp);

	free_aux_data(cd);
	free_diag_rho_data(cd);
}

int init_aux_data(ConfigData &cd)
{
	int num_global_states = pow(2, cd.Nc);

	cd.adjacement	= new int[num_global_states];
	cd.x_to_id		= new int[num_global_states];
	cd.id_to_x		= new int[cd.Ns];

	int state_id = 0;
	for (int global_state_id = 0; global_state_id < num_global_states; global_state_id ++)
	{
		if ((bit_count(global_state_id) == 2) && (bit_count(global_state_id & (global_state_id << 1)) == 1))
		{
			cd.adjacement[global_state_id] = 1;
		}
		else
		{
			cd.adjacement[global_state_id] = 0;
		}

		if (bit_count(global_state_id) == cd.Np)
		{
			cd.x_to_id[global_state_id] = state_id + 1;
			cd.id_to_x[state_id] = global_state_id;
			state_id ++;
		}
		else
		{
			cd.x_to_id[global_state_id] = 0;
		}
	}

	//print_int_array(cd.adjacement, num_global_states);
	//print_int_array(cd.x_to_id, num_global_states);
	//print_int_array(cd.id_to_x, cd.Ns);

	return state_id;
}

void free_aux_data(ConfigData &cd)
{
	delete[] cd.adjacement;
	delete[] cd.x_to_id;
	delete[] cd.id_to_x;
}

void init_hamiltonian_data(ConfigData &cd, ConfigParam &cp)
{
	double time = omp_get_wtime();

	cd.hamiltonian = new double[cd.Ns * cd.Ns];

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1 ++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2 ++)
		{
			cd.hamiltonian[state_id_1 * cd.Ns + state_id_2] = 0.0;
		}
	}

	init_hamiltonian_disorder_data(cd, cp);
	init_hamiltonian_interaction_data(cd, cp);
	init_hamiltonian_hopping_data(cd, cp);

	if (cp.dump >= 2)
	{
		string hamiltonian_fn = cp.path + "hamiltonian" + file_name_suffix(cp, 4);
		cout << "save hamiltonian to file:" << endl << hamiltonian_fn << endl << endl;
		write_double_data(hamiltonian_fn, cd.hamiltonian, cd.Ns * cd.Ns, 16);
	}

	time = omp_get_wtime() - time;
	cout << "time of init hamiltonian: " << time << endl << endl;
}

void init_hamiltonian_disorder_data(ConfigData &cd, ConfigParam &cp)
{
	double * energies = new double[cd.Nc];
	ifstream energies_ifs(cp.init);
	if (energies_ifs.good())
	{
		cout << "reading energies from file: " << endl << cp.init << endl << endl;

		int read_lines_count = 0;
		while (energies_ifs.good())
		{
			if (read_lines_count > cd.Nc)
			{
				stringstream msg;
				msg << "wrong file with energies: num lines in file more than num cells " << cd.Nc << endl;
				energies_ifs.close();
				Error(msg.str());
			}

			energies_ifs >> energies[read_lines_count];

			read_lines_count ++;
		}

		if (read_lines_count < cd.Nc)
		{
			stringstream msg;
			msg << "wrong file with energies: num lines in file less than num cells " << cd.Nc << endl;
			energies_ifs.close();
			Error(msg.str());
		}

		energies_ifs.close();
	}
	else
	{
		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, cp.seed, cp.max_num_seeds);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, cd.Nc, energies, -1.0, 1.0);

		string energies_fn = cp.path + "energies" + file_name_suffix(cp, 4);
		cout << "generate energies" << endl << "save energies to file:" << endl << energies_fn << endl << endl;
		write_double_data(energies_fn, energies, cd.Nc, 16);
	}

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		vector<int> vb = convert_int_to_vector_of_bits(cd.id_to_x[state_id], cd.Nc);
		double sum = 0.0;
		for (int cell_id = 0; cell_id < cd.Nc; cell_id ++)
		{
			sum += double(vb[cell_id]) * energies[cell_id];
		}
		sum *= 2.0 * cp.W;

		cd.hamiltonian[state_id * cd.Ns + state_id] += sum;
	}

	delete[] energies;
}

void init_hamiltonian_interaction_data(ConfigData &cd, ConfigParam &cp)
{
	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		cd.hamiltonian[state_id * cd.Ns + state_id] += cp.U * bit_count(cd.id_to_x[state_id] & (cd.id_to_x[state_id] << 1));
	}
}

void init_hamiltonian_hopping_data(ConfigData &cd, ConfigParam &cp)
{
	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
		{
			cd.hamiltonian[state_id_1 * cd.Ns + state_id_2] -= cp.J * cd.adjacement[cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2]];	
		}
	}
}

void free_hamiltonian_data(ConfigData &cd)
{
	delete[] cd.hamiltonian;
}

void init_hamiltonian_eig_data(ConfigData &cd, ConfigParam &cp)
{
	double time = omp_get_wtime();

	cd.hamiltonian_ev	= new double[cd.Ns * cd.Ns];
	cd.hamiltonian_eg	= new double[cd.Ns];

	int founded_eigen_vals = 0;
	int * isuppz = new int[2 * cd.Ns];

	int info = LAPACKE_dsyevr(
		LAPACK_ROW_MAJOR,
		'V',
		'A',
		'U',
		cd.Ns,
		cd.hamiltonian,
		cd.Ns,
		0,
		0,
		0,
		0,
		EPS,
		&founded_eigen_vals,
		cd.hamiltonian_eg,
		cd.hamiltonian_ev,
		cd.Ns,
		isuppz
		);

	delete[] isuppz;

	if( info > 0 ) 
	{
		stringstream msg;
		msg << "failed to solve hamiltonian eigenproblem "<< endl;
		Error(msg.str());
	}

	if (cp.dump >= 1)
	{
		string hamiltonian_eg_fn = cp.path + "hamiltonian_eg" + file_name_suffix(cp, 4);
		cout << "save hamiltonian eigen values to file:" << endl << hamiltonian_eg_fn << endl << endl;
		write_double_data(hamiltonian_eg_fn, cd.hamiltonian_eg, cd.Ns, 16);
	}

	if (cp.dump >= 2)
	{
		string hamiltonian_ev_fn = cp.path + "hamiltonian_ev" + file_name_suffix(cp, 4);
		cout << "save hamiltonian eigen vectors to file:" << endl << hamiltonian_ev_fn << endl << endl;
		write_double_data(hamiltonian_ev_fn, cd.hamiltonian_ev, cd.Ns * cd.Ns, 16);
	}

	time = omp_get_wtime() - time;
	cout << "time of solving hamiltonian eigenproblem: " << time << endl << endl;
}

void free_hamiltonian_eig_data(ConfigData &cd)
{
	delete[] cd.hamiltonian_ev;
	delete[] cd.hamiltonian_eg;
}

void init_trans_rates_data(ConfigData &cd, ConfigParam &cp)
{
	double time = omp_get_wtime();

	cd.trans_rates = new double[cd.Ns * cd.Ns];
	
	MKL_Complex16 * hamiltonian_ev_complex = new MKL_Complex16[cd.Ns * cd.Ns];

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
		{
			cd.trans_rates[state_id_1 * cd.Ns + state_id_2] = 0.0;

			hamiltonian_ev_complex[state_id_1 * cd.Ns + state_id_2].real = cd.hamiltonian_ev[state_id_1 * cd.Ns + state_id_2];
			hamiltonian_ev_complex[state_id_1 * cd.Ns + state_id_2].imag = 0.0;
		}
	}

	MKL_Complex16 * tmp_1 = new MKL_Complex16[cd.Ns * cd.Ns]; // dissipator and final result 
	MKL_Complex16 * tmp_2 = new MKL_Complex16[cd.Ns * cd.Ns]; // intermediate result

	for (int dissipator_id = 0; dissipator_id < cd.Nc - 1; dissipator_id ++)
	{
		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
			{
				tmp_1[state_id_1 * cd.Ns + state_id_2].real = 0.0;
				tmp_1[state_id_1 * cd.Ns + state_id_2].imag = 0.0;

				tmp_2[state_id_1 * cd.Ns + state_id_2].real = 0.0;
				tmp_2[state_id_1 * cd.Ns + state_id_2].imag = 0.0;
			}
		}

		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			tmp_1[state_id_1 * cd.Ns + state_id_1].real = double(bit_at(cd.id_to_x[state_id_1], dissipator_id)) - double(bit_at(cd.id_to_x[state_id_1], dissipator_id + 1));
			tmp_1[state_id_1 * cd.Ns + state_id_1].imag = 0.0;

			for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
			{
				if (cd.adjacement[cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2]])
				{
					vector<int> adjacency_bits = convert_int_to_vector_of_bits(cd.id_to_x[state_id_1] ^ cd.id_to_x[state_id_2], cd.Nc);
					vector<int> hop;
					for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
					{
						if(adjacency_bits[cell_id])
						{
							hop.push_back(cell_id);
						}
					}

					for (int ad_cell_id = 0; ad_cell_id < hop.size(); ad_cell_id++)
					{
						hop[ad_cell_id] = (cd.Nc - 1) - hop[ad_cell_id];
					}

					if (hop[1] == dissipator_id)
					{
						if (bit_at(cd.id_to_x[state_id_1], dissipator_id))
						{
							tmp_1[state_id_2 * cd.Ns + state_id_1].real = cos(cp.alpha);
							tmp_1[state_id_2 * cd.Ns + state_id_1].imag = sin(cp.alpha);

							tmp_1[state_id_1 * cd.Ns + state_id_2].real = -cos(-cp.alpha);
							tmp_1[state_id_1 * cd.Ns + state_id_2].imag = -sin(-cp.alpha);
						}
						else
						{
							tmp_1[state_id_2 * cd.Ns + state_id_1].real = -cos(-cp.alpha);
							tmp_1[state_id_2 * cd.Ns + state_id_1].imag = -sin(-cp.alpha);

							tmp_1[state_id_1 * cd.Ns + state_id_2].real = cos(cp.alpha);
							tmp_1[state_id_1 * cd.Ns + state_id_2].imag = sin(cp.alpha);
						}
					}

					adjacency_bits.clear();
					hop.clear();
				}
			}
		}

		MKL_Complex16 ZERO	= {0.0, 0.0};
		MKL_Complex16 ONE	= {1.0, 0.0};

		cblas_zgemm(
			CblasRowMajor,
			CblasConjTrans,
			CblasNoTrans,
			cd.Ns,
			cd.Ns,
			cd.Ns,
			&ONE,
			hamiltonian_ev_complex,
			cd.Ns,
			tmp_1,
			cd.Ns,
			&ZERO,
			tmp_2,
			cd.Ns
			);

		cblas_zgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			cd.Ns,
			cd.Ns,
			cd.Ns,
			&ONE,
			tmp_2,
			cd.Ns,
			hamiltonian_ev_complex,
			cd.Ns,
			&ZERO,
			tmp_1,
			cd.Ns
			);

		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
			{
				cd.trans_rates[state_id_1 * cd.Ns + state_id_2] += cp.g * (tmp_1[state_id_1 * cd.Ns + state_id_2].real * tmp_1[state_id_1 * cd.Ns + state_id_2].real + tmp_1[state_id_1 * cd.Ns + state_id_2].imag * tmp_1[state_id_1 * cd.Ns + state_id_2].imag);
			}
		}
	}

	delete[] hamiltonian_ev_complex;

	delete[] tmp_2;
	delete[] tmp_1;


	double * column_sums = new double[cd.Ns];
	for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
	{
		column_sums[state_id_2] = 0.0;

		for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
		{
			column_sums[state_id_2] += cd.trans_rates[state_id_1 * cd.Ns + state_id_2];
		}

		cd.trans_rates[state_id_2 * cd.Ns + state_id_2] -= column_sums[state_id_2];
	}
	delete[] column_sums;

	if (cp.dump >= 2)
	{
		string trans_rates_fn = cp.path + "trans_rates" + file_name_suffix(cp, 4);
		cout << "save trans rates to file:" << endl << trans_rates_fn << endl << endl;
		write_double_data(trans_rates_fn, cd.trans_rates, cd.Ns * cd.Ns, 16);
	}

	time = omp_get_wtime() - time;
	cout << "time of creating trans rates matrix: " << time << endl << endl;
}

void free_trans_rates_data(ConfigData &cd)
{
	delete[] cd.trans_rates;
}

void init_diag_rho_data(ConfigData &cd, ConfigParam &cp)
{
	double time = omp_get_wtime();

	cd.diag_rho_in_st = new double[cd.Ns];
	cd.diag_rho_in_d = new double[cd.Ns];
	cd.rho_in_d = new double[cd.Ns * cd.Ns];

	double * eg_real = new double[cd.Ns];
	double * eg_imag = new double[cd.Ns];
	double * ev = new double[cd.Ns * cd.Ns];

	int info = LAPACKE_dgeev( 
		LAPACK_ROW_MAJOR,
		'N',
		'V',
		cd.Ns,
		cd.trans_rates,
		cd.Ns,
		eg_real,
		eg_imag,
		NULL,
		cd.Ns,
		ev,
		cd.Ns
		);

	if( info > 0 ) 
	{
		stringstream msg;
		msg << "failed to solve hamiltonian eigenproblem "<< endl;
		Error(msg.str());
	}

	vector<double> eg_reals_vec;
	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		eg_reals_vec.push_back(fabs(eg_real[state_id]));
	}

	vector<int> order = sort_doubles_with_order(eg_reals_vec);

	eg_reals_vec.clear();

	if (eg_imag[order[0]] > EPS)
	{
		cout << "zero eg is complex: " << eg_real[order[0]] << " + " << eg_imag[order[0]] << " i" << endl;
		stringstream msg;
		msg << "diag rho must be non-complex "<< endl;
		Error(msg.str());
	}


	double sum = 0.0;
	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		cd.diag_rho_in_st[state_id] = ev[state_id * cd.Ns + order[0]];
		sum += cd.diag_rho_in_st[state_id];
	}

	order.clear();

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		cd.diag_rho_in_st[state_id] /= sum;
	}

	if (cp.dump >= 1)
	{
		string diag_rho_in_st_fn = cp.path + "diag_rho_in_st" + file_name_suffix(cp, 4);
		cout << "save diag rho in stationary basis to file:" << endl << diag_rho_in_st_fn << endl << endl;
		write_double_data(diag_rho_in_st_fn, cd.diag_rho_in_st, cd.Ns, 16);
	}

	delete[] ev;
	delete[] eg_real;
	delete[] eg_imag;

	double * tmp_1 = new double[cd.Ns * cd.Ns]; // rho diag in stationary basis and final result 
	double * tmp_2 = new double[cd.Ns * cd.Ns]; // intermediate result

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
		{
			tmp_1[state_id_1 * cd.Ns + state_id_2] = 0.0;
			tmp_2[state_id_1 * cd.Ns + state_id_2] = 0.0;
		}
	}

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		tmp_1[state_id * cd.Ns + state_id] = cd.diag_rho_in_st[state_id];
	}

	double ZERO	= 0.0;
	double ONE = 1.0;

	cblas_dgemm(
		CblasRowMajor,
		CblasNoTrans,
		CblasNoTrans,
		cd.Ns,
		cd.Ns,
		cd.Ns,
		ONE,
		cd.hamiltonian_ev,
		cd.Ns,
		tmp_1,
		cd.Ns,
		ZERO,
		tmp_2,
		cd.Ns
		);

	cblas_dgemm(
		CblasRowMajor,
		CblasNoTrans,
		CblasConjTrans,
		cd.Ns,
		cd.Ns,
		cd.Ns,
		ONE,
		tmp_2,
		cd.Ns,
		cd.hamiltonian_ev,
		cd.Ns,
		ZERO,
		tmp_1,
		cd.Ns
		);

	for (int state_id_1 = 0; state_id_1 < cd.Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < cd.Ns; state_id_2++)
		{
			cd.rho_in_d[state_id_1 * cd.Ns + state_id_2] = tmp_1[state_id_1 * cd.Ns + state_id_2];
		}

		cd.diag_rho_in_d[state_id_1] = tmp_1[state_id_1 * cd.Ns + state_id_1];	
	}

	delete[] tmp_1;
	delete[] tmp_2;

	if (cp.dump >= 1)
	{
		string diag_rho_in_d_fn = cp.path + "diag_rho_in_d" + file_name_suffix(cp, 4);
		cout << "save diag rho in direct basis to file:" << endl << diag_rho_in_d_fn << endl << endl;
		write_double_data(diag_rho_in_d_fn, cd.diag_rho_in_d, cd.Ns, 16);
	}

	if (cp.dump >= 2)
	{
		string rho_in_d_fn = cp.path + "rho_in_d" + file_name_suffix(cp, 4);
		cout << "save rho in direct basis to file:" << endl << rho_in_d_fn << endl << endl;
		write_double_data(rho_in_d_fn, cd.rho_in_d, cd.Ns * cd.Ns, 16);
	}

	time = omp_get_wtime() - time;
	cout << "time of calculating rho: " << time << endl << endl;
}

void free_diag_rho_data(ConfigData &cd)
{
	delete[] cd.diag_rho_in_st;
	delete[] cd.rho_in_d;
	delete[] cd.diag_rho_in_d;
}

void calculate_characteristics(ConfigData &cd, ConfigParam &cp)
{
	double time = omp_get_wtime();

	// ######## entropy ########
	cd.entropy = 0.0;
	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		cd.entropy -= cd.diag_rho_in_st[state_id] * log(cd.diag_rho_in_st[state_id]);
	}

	// ######## imbalance ########
	double * n_part = new double[cd.Nc];

	for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
	{
		n_part[cell_id] = 0.0;
	}

	for (int state_id = 0; state_id < cd.Ns; state_id++)
	{
		vector<int> vb = convert_int_to_vector_of_bits(cd.id_to_x[state_id], cd.Nc);

		for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
		{
			n_part[cell_id] += cd.rho_in_d[state_id * cd.Ns + state_id] * double(vb[cell_id]);
		}
	}

	double sum_odd = 0.0;
	double sum_even = 0.0;
	double sum_all = 0.0;


	for (int cell_id = 0; cell_id < cd.Nc; cell_id++)
	{
		if (cell_id%2 == 0)
		{
			sum_odd += n_part[cell_id];
		}
		else
		{
			sum_even += n_part[cell_id];
		}
	}
	sum_all = sum_even + sum_odd;

	cd.imbalance = (sum_odd - sum_even) / sum_all;

	delete[] n_part;

	double * characteristics = new double[2];
	characteristics[0] = cd.entropy;
	characteristics[1] = cd.imbalance;

	string characteristics_fn = cp.path + "characteristics" + file_name_suffix(cp, 4);
	cout << "save characteristics to file:" << endl << characteristics_fn << endl << endl;
	write_double_data(characteristics_fn, characteristics, 2, 16);

	delete[] characteristics;

	time = omp_get_wtime() - time;
	cout << "time of calculating characteristics: " << time << endl << endl;
}