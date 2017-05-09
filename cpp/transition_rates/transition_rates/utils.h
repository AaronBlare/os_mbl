#pragma once
#include "config.h"

#define Error( msg ) error( msg, __FUNCTION__, __FILE__, __LINE__ )

#define PI 3.1415926535897932384626433832795

#define EPS 1.0e-12

void error(const std::string& err, const char* func, const char* file, int line);

int n_choose_k(int n, int k);

int bit_count(int value);

int bit_at(int value, int position);

void print_int_array(int * data, int N);

string file_name_suffix(ConfigParam &cp, int precision);

void write_double_data(string file_name, double * data, int size, int precision);

void write_complex_data(string file_name, MKL_Complex16 * data, int size, int precision);

vector<int> convert_int_to_vector_of_bits(int x, int size);

vector<int> sort_doubles_with_order(vector<double> &v);
