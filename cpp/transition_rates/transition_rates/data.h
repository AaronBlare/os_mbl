#pragma once
#include "config.h"
#include "utils.h"


struct ConfigData
{
	int Nc;		// number of cells
	int Np;		// number of particles
	int Ns;		// number of states

	int * adjacement;	// aux array for define adjacement states
	int * x_to_id;		// aux array for id definition from binary representation
	int * id_to_x;		// aux array for binary representation definition from id

	double * hamiltonian;		// hamiltonian
	double * hamiltonian_ev;		// eigen vectors of hamiltonian
	double * hamiltonian_eg;		// eigen values of hamiltonian

	double * trans_rates;	// transition rates matrix

	MKL_Complex16 * lindbladian;		// libladian

	double * diag_rho_in_st;			// diag of matrix rho in stationary basis
	double * rho_in_d;			// rho in direct basis
	double * diag_rho_in_d;				// diag of matrix rho in direct basis
	double entropy;						// entropy of asymptotic state
	double imbalance;					// entropy of asymptotic state
};

void run_trans_rates(ConfigData &cd, ConfigParam &cp);
void run_zero_eigen_vector(ConfigData &cd, ConfigParam &cp);

int init_aux_data(ConfigData &cd);
void free_aux_data(ConfigData &cd);

void init_hamiltonian_data(ConfigData &cd, ConfigParam &cp);
void init_hamiltonian_disorder_data(ConfigData &cd, ConfigParam &cp);
void init_hamiltonian_interaction_data(ConfigData &cd, ConfigParam &cp);
void init_hamiltonian_hopping_data(ConfigData &cd, ConfigParam &cp);
void free_hamiltonian_data(ConfigData &cd);

void init_hamiltonian_eig_data(ConfigData &cd, ConfigParam &cp);
void free_hamiltonian_eig_data(ConfigData &cd);

void init_trans_rates_data(ConfigData &cd, ConfigParam &cp);
void free_trans_rates_data(ConfigData &cd);

void init_diag_rho_data(ConfigData &cd, ConfigParam &cp);
void free_diag_rho_data(ConfigData &cd);

void calculate_characteristics(ConfigData &cd, ConfigParam &cp);

void init_libladian(ConfigData &cd, ConfigParam &cp);
void free_libladian(ConfigData &cd, ConfigParam &cp);
