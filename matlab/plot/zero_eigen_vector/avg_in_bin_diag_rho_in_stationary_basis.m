clear all;

W = 8; 	% disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

Nc = 8;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

diss_type = 1; 		% Dissipator type: 0-Poletti, 1-Diehl
alpha = pi; 			% Dissipator Diehl phase 
energy_type = 0;   	% 0 if regular, 1 if zero mean

start_seed_start = 1;
shift_seed_start = 100;
num_seed_start = 100;

num_seeds = shift_seed_start * num_seed_start;

num_intervals_eg = 100;

size_total = num_seeds * Ns;

egs         = zeros(size_total, 1);
diag_rho    = zeros(size_total, 1);

for seed_start = start_seed_start : shift_seed_start: start_seed_start + shift_seed_start * (num_seed_start - 1)

    seed_start = seed_start
  
    path = sprintf('../results/dt_%d/alpha_%0.4f/et_%d/Ns_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/seed_start_%d', ...
        diss_type,alpha, energy_type, Nc, W, U, J, g, seed_start);
    
    for seed = seed_start : seed_start + (shift_seed_start-1)
	
        file_name = sprintf('%s/eg_hamiltonian_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', ...
            path, Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
			
        local_data_eg = importdata(file_name);
		
		file_name = sprintf('%s/rho_diag_in_stationary_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', ...
            path, Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
		
		local_data_rho = importdata(file_name);
        
	    for state_id = 1:Ns
			egs((seed - 1) * Ns + state_id) = local_data_eg(state_id);
			diag_rho((seed - 1) * Ns + state_id) = local_data_rho(state_id);
		end
		
    end
end

min_eg = min(egs)
max_eg = max(egs)

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

savefig(sprintf('avg_in_bin_diag_rho_in_stationary_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ss(%d)_sn(%d).fig', ...
	Nc, diss_type, alpha, energy_type, W, U, J, g, seed_start, num_seeds));