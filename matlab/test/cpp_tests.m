clear all;

filename = '../source/int/config.txt';
input_data = importdata(filename);

Nc                  = input_data(1);
diss_type           = input_data(2);
diss_phase          = input_data(3);
energy_type         = input_data(4);
periodic_bc 		= input_data(5);
W                   = input_data(6);
U                   = input_data(7);
J                   = input_data(8);
g                   = input_data(9);
seed_start          = input_data(10);
seed_num            = input_data(11);
is_int              = input_data(12);
int_ist             = input_data(13);
int_isi             = input_data(14);
int_dt              = input_data(15);
int_db              = input_data(16);
int_de              = input_data(17);
int_dn              = input_data(18);
is_zev              = input_data(19);
is_zev_check        = input_data(20);
is_eg_stat          = input_data(21);
is_save_vec         = input_data(22);
is_save_mtx         = input_data(23);
fs_type             = input_data(24);


data_path = '../../data/matlab/';

Np = Nc/2;
Ns = nchoosek(Nc, Np);

seed = seed_start;

fn_suffix = sprintf('zev_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    periodic_bc, ...
    W, ...
    U, ...
    J, ...
    g, ...
    seed);

file_name = sprintf('%srho_d_%s', data_path, fn_suffix);
rho_d_data = importdata(file_name);
rho_d_zev = zeros(Ns, Ns);
for s_id_1 = 1:Ns
    for s_id_2 = 1:Ns
        rho_d_zev(s_id_1, s_id_2) = rho_d_data((s_id_1 - 1) * Ns + s_id_2, 1) + sqrt(-1) * rho_d_data((s_id_1 - 1) * Ns + s_id_2, 2);
    end
end

data_path = '../../cpp/os_mbl/os_mbl/';
fn_suffix = sprintf('f_int_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(0)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(1).txt', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    W, ...
    U, ...
    J, ...
    g);

file_name = sprintf('%srho_in_d_%s', data_path, fn_suffix);
rho_d_data = importdata(file_name);

rho_d_fb = zeros(Ns, Ns);

for i = 1:Ns
    for j = 1:Ns
        rho_d_fb(i,j) = rho_d_data((i-1) * Ns + j, 1) + sqrt(-1) * rho_d_data((i-1) * Ns + j, 2);
    end
end

diff_mtx = rho_d_fb - rho_d_zev;

max_diff = max(max(abs(diff_mtx)))


Nss = 2^Nc;

file_name = sprintf('%srho_xtd_%s', data_path, fn_suffix);
rho_xtd_data = importdata(file_name);

rho_xtd_fb = zeros(Nss, Nss);

for i = 1:Nss
    for j = 1:Nss
        rho_xtd_fb(i,j) = rho_xtd_data((i-1) * Nss + j, 1) + sqrt(-1) * rho_xtd_data((i-1) * Nss + j, 2);
    end
end

redNss = 2^(Nc/2);

file_name = sprintf('%srho_red_%s', data_path, fn_suffix);
rho_red_data = importdata(file_name);

rho_red_fb = zeros(redNss, redNss);

for i = 1:redNss
    for j = 1:redNss
        rho_red_fb(i,j) = rho_red_data((i-1) * redNss + j, 1) + sqrt(-1) * rho_red_data((i-1) * redNss + j, 2);
    end
end



