clear all;

Nc = 14;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

dt = 1; 		% Dissipator type: 0-Poletti, 1-Diehl
alpha = pi; 			% Dissipator Diehl phase 
et = 0;   	% 0 if regular, 1 if zero mean

W = 8; 	% disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

max_num_seeds = 1000000;
seed_start = 1;
num_seeds = 200;

num_intervals_eg = 100;
size_del_tails = 0;

size_total = num_seeds * Ns;

egs         = zeros(size_total, 1);
diag_rho    = zeros(size_total, 1);

for seed = seed_start : seed_start + (num_seeds - 1)

	seed = seed
	start_id = (seed - seed_start) * Ns + 1
	
    path = sprintf('../results/Nc_%d/dt_%d/alpha_%0.4f/et_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/max_num_seeds_%d/seed_%d', ...
        Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
    
    file_name = sprintf('%s/hamiltonian_eg_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', ...
        path, Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
		
	file_name = file_name
    data = importdata(file_name);
    for state_id = 1:Ns
        egs((seed - seed_start) * Ns + state_id) = data(state_id);
    end

    file_name = sprintf('%s/diag_rho_in_st_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', ...
        path, Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed); 
    data = importdata(file_name);
    for state_id = 1:Ns
        diag_rho((seed - seed_start) * Ns + state_id) = data(state_id);
    end
end

min_eg = min(egs);
max_eg = max(egs);

shift_eg = (max_eg - min_eg) / num_intervals_eg;
intervals_eg = zeros(num_intervals_eg, 1);
diag_rho_in_intervals_avg = zeros(num_intervals_eg, 1);
count_in_intervals = zeros(num_intervals_eg, 1);
for i = 1:num_intervals_eg
    intervals_eg(i) = min_eg + i * shift_eg - 0.5 * shift_eg;
end

for i = 1:size_total
    
    current_eg = egs(i);
    current_diag_rho = diag_rho(i); 
    
    if ((current_eg >= min_eg) && (current_eg <= max_eg))
        id_eg = floor((current_eg - min_eg) * num_intervals_eg / (max_eg - min_eg + 0.000001)) + 1;
        diag_rho_in_intervals_avg(id_eg) = diag_rho_in_intervals_avg(id_eg) + diag_rho(i);
        count_in_intervals(id_eg) = count_in_intervals(id_eg) + 1;
    end
end

del_indexes = horzcat([1:size_del_tails], [num_intervals_eg - size_del_tails + 1: num_intervals_eg])

diag_rho_in_intervals_avg(del_indexes) = [];
count_in_intervals(del_indexes) = [];
intervals_eg(del_indexes) = [];

num_intervals_eg = num_intervals_eg

num_intervals_eg = num_intervals_eg - 2 * size_del_tails;

num_intervals_eg = num_intervals_eg

min_eg = min(intervals_eg);
max_eg = max(intervals_eg);

for i = 1:num_intervals_eg
    diag_rho_in_intervals_avg(i) = diag_rho_in_intervals_avg(i) / count_in_intervals(i);
	intervals_eg(i) = (intervals_eg(i) - min_eg) / (max_eg - min_eg) * 2 - 1;
end

norm = 0.0;
for i = 1:num_intervals_eg
    norm = norm + diag_rho_in_intervals_avg(i);
end

for i = 1:num_intervals_eg
    diag_rho_in_intervals_avg(i) = diag_rho_in_intervals_avg(i) / norm;
end

figure;
hLine = plot(intervals_eg, diag_rho_in_intervals_avg, '-o', 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$\epsilon$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$<|\rho_{p,p}|>$', 'Interpreter', 'latex');

savefig(sprintf('avg_in_bin_diag_rho_in_stationary_basis_normalized_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_ss(%d)_sn(%d).fig', ...
	Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed_start, num_seeds));
