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

states_num = nchoosek(Nc, Nc/2);

type = 2;
begin_eval_id = floor(states_num * 1/3);
end_eval_id = floor(states_num * 2/3);
begin_part = 1 / 3;
end_part = 2 / 3;
begin_eval = 1 * 10.^(-6);
end_eval = 1 * 10.^(-4);

num_int = 1000;
evals_diff_begin = 0.0;
evals_diff_end = 10.0;
evals_diff_shift = (evals_diff_end - evals_diff_begin) / num_int;
evals_diff_int = zeros(num_int, 1);
evals_diff_pdf = zeros(num_int, 1);
for int_id = 1:num_int
    evals_diff_int(int_id) = evals_diff_begin + int_id * evals_diff_shift - 0.5 * evals_diff_shift;
end
eps = 1.0e-8;

total_evals_num = 0;
total_evals_diff_num = 0;

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
        
        file_name = sprintf('%s/rho_evals_%s', curr_path, suffix);
        
        curr_evals = importdata(file_name);
        curr_evals = sort(curr_evals);
        
        
        if type == 0
            
            curr_evals = curr_evals(begin_eval_id:end_eval_id);
            num_evals = size(curr_evals, 1)
            
            mean_diff = (curr_evals(num_evals) - curr_evals(1)) / num_evals
            
            for st_id = 1:num_evals - 1
                curr_eval_diff = (curr_evals(st_id + 1) - curr_evals(st_id)) / mean_diff;
                
                if (curr_eval_diff < evals_diff_end)
                    int_id = floor((curr_eval_diff - evals_diff_begin) * num_int / (evals_diff_end - evals_diff_begin + eps)) + 1;
                    evals_diff_pdf(int_id) = evals_diff_pdf(int_id) + 1;
                end
            end
            
        elseif type == 1
            
            curr_evals(1)
            total_int = curr_evals(end) - curr_evals(1)
            begin_lim = curr_evals(1) + total_int * begin_part
            end_lim = curr_evals(1) + total_int * end_part
            
            curr_eval_id = 1;
            
            while curr_evals(curr_eval_id) < begin_lim
                curr_eval_id = curr_eval_id + 1;
            end
            
            begin_eval_id = curr_eval_id
            
            while curr_evals(curr_eval_id) < end_lim
                curr_eval_id = curr_eval_id + 1;
            end
            
            end_eval_id = curr_eval_id - 1
            
            curr_evals = curr_evals(begin_eval_id:end_eval_id);
            num_evals = size(curr_evals, 1);
            
            if (num_evals >= 2)
			
				total_evals_diff_num = total_evals_diff_num + (num_evals - 1);
                
                mean_diff = (curr_evals(num_evals) - curr_evals(1)) / num_evals;
                
                for st_id = 1:num_evals - 1
                    curr_eval_diff = (curr_evals(st_id + 1) - curr_evals(st_id)) / mean_diff;
                    
                    if (curr_eval_diff < evals_diff_end)
                        int_id = floor((curr_eval_diff - evals_diff_begin) * num_int / (evals_diff_end - evals_diff_begin + eps)) + 1;
                        evals_diff_pdf(int_id) = evals_diff_pdf(int_id) + 1;
                    end
                end
                
            end
            
        elseif type == 2
            
            begin_lim = begin_eval;
            end_lim = end_eval;
            
            curr_eval_id = 1;
            
            while curr_evals(curr_eval_id) < begin_lim
                curr_eval_id = curr_eval_id + 1;
            end
            
            begin_eval_id = curr_eval_id;
            
            while curr_evals(curr_eval_id) < end_lim
                curr_eval_id = curr_eval_id + 1;
            end
            
            end_eval_id = curr_eval_id - 1;
            
            curr_evals = curr_evals(begin_eval_id:end_eval_id);
            num_evals = size(curr_evals, 1);
			
			total_evals_num = total_evals_num + num_evals;
            
            if (num_evals >= 2)
			
				total_evals_diff_num = total_evals_diff_num + (num_evals - 1);
                
                mean_diff = (curr_evals(num_evals) - curr_evals(1)) / num_evals;
                
                for st_id = 1:num_evals - 1
                    curr_eval_diff = (curr_evals(st_id + 1) - curr_evals(st_id)) / mean_diff;
                    
                    if (curr_eval_diff < evals_diff_end)
                        int_id = floor((curr_eval_diff - evals_diff_begin) * num_int / (evals_diff_end - evals_diff_begin + eps)) + 1;
                        evals_diff_pdf(int_id) = evals_diff_pdf(int_id) + 1;
                    end
                end
                
            end
            
        end
        
    end
end

if type == 0
    evals_diff_pdf = evals_diff_pdf / (seed_num * seed_start_num * (num_evals-1) * evals_diff_shift);
end

if type > 0
    evals_diff_pdf = evals_diff_pdf / (evals_diff_shift * total_evals_diff_num);
end

total_evals_num = total_evals_num

sum(evals_diff_pdf)

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

file_name = sprintf('rho_evals_diff_pdf_%s', fn_suffix);
file_id = fopen(file_name, 'w');
for int_id = 1:num_int
    fprintf(file_id, '%0.18e %0.18e\n', evals_diff_int(int_id), evals_diff_pdf(int_id));
end
fclose(file_id);


