clear all;

Nc = 8;
W = 8.0;
U = 1.0;
J = 1.0;
g = 0.1;
dt = 1;
alpha = 1.0;
et = 0;
seed_matlab = 1337;
seed_cpp = 1;
max_num_seeds = 1000;

Np = Nc/2;
Ns = nchoosek(Nc,Np); 

path_matlab = '../../data/transition_rates/matlab/';
path_cpp = '../../data/transition_rates/cpp/';

suffix_matlab = sprintf('_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, dt, alpha, et, W, U, J, g, seed_matlab);
suffix_cpp = sprintf('_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed_cpp);

file_name = sprintf('%shamiltonian%s', path_matlab, suffix_matlab);
data = importdata(file_name);
mtx_matlab = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

file_name = sprintf('%shamiltonian%s', path_cpp, suffix_cpp);
data = importdata(file_name);
mtx_cpp = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

H_diff = max(max(abs(mtx_matlab - mtx_cpp)))


file_name = sprintf('%stransition_rates%s', path_matlab, suffix_matlab);
data = importdata(file_name);
mtx_matlab = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

file_name = sprintf('%strans_rates%s', path_cpp, suffix_cpp);
data = importdata(file_name);
mtx_cpp = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

transition_rates_diff = max(max(abs(mtx_matlab - mtx_cpp)))
