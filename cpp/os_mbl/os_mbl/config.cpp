#include "config.h"

vector<string> split(const string& str, const string& delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);

		if (pos == string::npos) 
		{
			pos = str.length();
		}

		string token = str.substr(prev, pos-prev);

		if (!token.empty())
		{
			tokens.push_back(token);
		}

		prev = pos + delim.length();
	}
	while (pos < str.length() && prev < str.length());

	return tokens;
}

void set_param(ConfigParam &param, string str, string value)
{
	if (str.compare("task") == 0)
	{
		param.task = atoi(value.c_str());
	}

	if (str.compare("Nc") == 0)
	{
		param.Nc = atoi(value.c_str());
	}

	if (str.compare("dt") == 0)
	{
		param.dt = atoi(value.c_str());
	}

	if (str.compare("dp") == 0)
	{
		param.dp = atof(value.c_str());
	}

	if (str.compare("et") == 0)
	{
		param.et = atoi(value.c_str());
	}

	if (str.compare("bc") == 0)
	{
		param.bc = atoi(value.c_str());
	}

	if (str.compare("W") == 0)
	{
		param.W = atof(value.c_str());
	}

	if (str.compare("U") == 0)
	{
		param.U = atof(value.c_str());
	}
	
	if (str.compare("J") == 0)
	{
		param.J = atof(value.c_str());
	}

	if (str.compare("g") == 0)
	{
		param.g = atof(value.c_str());
	}

	if (str.compare("seed") == 0)
	{
		param.seed = atoi(value.c_str());
	}

	if (str.compare("max_num_seeds") == 0)
	{
		param.max_num_seeds = atoi(value.c_str());
	}

	if (str.compare("dump_vecs") == 0)
	{
		param.dump_vecs = atoi(value.c_str());
	}

	if (str.compare("dump_mtxs") == 0)
	{
		param.dump_mtxs = atoi(value.c_str());
	}

	if (str.compare("dump_sops") == 0)
	{
		param.dump_sops = atoi(value.c_str());
	}

	if (str.compare("init_file") == 0)
	{
		param.init_file = value;
	}

	if (str.compare("path") == 0)
	{
		param.path = value;
	}
}

void init_config_param(ConfigParam &param, char * file_name)
{
	string line;
	ifstream config_file (file_name);
	if (config_file.is_open())
	{
		while (getline(config_file, line))
		{
			vector<string> tokens = split(line, " ");

			if (tokens.size() == 2)
			{
				set_param(param, tokens[0], tokens[1]);
			}
		}
		config_file.close();
	}
	else 
	{
		cout << "unable to open file" << endl;
		cout << "init with default params" << endl;
	}

	output_setting(param);
}

void output_setting(ConfigParam &param)
{
	cout << "############# parameters #############"					<< endl;
	cout << "task = " << param.task										<< endl;
	cout << "Nc = " << param.Nc											<< endl;
	cout << "dt = " << param.dt											<< endl;
	cout << "dp = " << param.dp											<< endl;
	cout << "et = " << param.et											<< endl;
	cout << "bc = " << param.bc											<< endl;
	cout << "W = " << param.W											<< endl;
	cout << "U = " << param.U											<< endl;
	cout << "J = " << param.J											<< endl;
	cout << "g = " << param.g											<< endl;
	cout << "seed = " << param.seed										<< endl;

	cout << "max_num_seeds = " << param.max_num_seeds					<< endl;
	cout << "dump_vecs = " << param.dump_vecs							<< endl;
	cout << "dump_mtxs = " << param.dump_mtxs							<< endl;
	cout << "dump_sops = " << param.dump_sops							<< endl;
	cout << "init_file = " << param.init_file							<< endl;
	cout << "path = " << param.path										<< endl;
	cout << "######################################"					<< endl;
}