clear all;


data_path = '/data/biophys/yusipov/os_mbl/int/matlab/';
prefix = 'characteristics';
data_path = sprintf('%s%s', data_path, prefix);

filename = 'config.txt';
input_data = importdata(filename);


Nc = 8;
diss_type = 1;
diss_phase = 0.0;
energy_type = 0;
periodic_bc = 0;
W = input_data(1);
U = 1.0;
J = 1.0;
g = 0.1;
seed_start = 1;
seed_num = 100;
is_int = 0;
int_ist = 0;
int_isi = 50;
int_dt = 1;
int_db = 0.1;
int_de = 10000;
int_dn = 250;
is_zev = 1;
is_zev_check = 1;
is_eg_stat = 1;
is_save_vec = 1;
is_save_mtx = 0;
fs_type = 0;

seed_start_begin = seed_start;
seed_start_num = 100;

num_int = 100;
int_begin = -0.5;
int_end = 0.5;
int_shift = (int_end - int_begin) / num_int;
eps = 1.0e-6;

imb_int = zeros(num_int, 1);
imb_evo = zeros(num_int, 1);
for int_id = 1:num_int
    imb_int(int_id) = int_begin + int_id * int_shift - 0.5 * int_shift;
end
    
for seed_start = seed_start_begin : seed_num : seed_start_begin + seed_num * (seed_start_num - 1)
    
    for seed = seed_start : seed_start + (seed_num-1)
        
        curr_path = sprintf('%s/Nc_%d/dt_%d/dp_%0.4f/et_%d/bc_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/seed_start%d', ...
            data_path, ...
            Nc, ...
            diss_type, ...
            diss_phase, ...
            energy_type, ...
            periodic_bc, ...
            W, ...
            U, ...
            J, ...
            g, ...
            seed_start);
        
        suffix = sprintf('zev_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
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
        
        file_name = sprintf('%s/imbalance_%s', curr_path, suffix);
        
        imb_dst = importdata(file_name);
        
        int_id = floor((imb_dst - int_begin) * num_int / (int_end - int_begin + eps)) + 1;
        imb_evo(int_id) = imb_evo(int_id) + 1;
        
    end
end

imb_evo = imb_evo / (seed_num * seed_start_num * int_shift);

fn_suffix = sprintf('Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(var).txt', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    periodic_bc, ...
    W, ...
    U, ...
    J, ...
    g);

file_name = sprintf('imbalance_pdf_%s', fn_suffix);
file_id = fopen(file_name, 'w');
for int_id = 1:num_int
    fprintf(file_id, '%0.18e %0.18e\n', imb_int(int_id), imb_evo(int_id));
end
fclose(file_id);


