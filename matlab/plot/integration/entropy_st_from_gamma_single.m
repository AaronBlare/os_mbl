clear all;

W = 10; % disorder
U = 1;  % interaction
J = 1;  % hopping
g = 1;  % gamma;

Nc = 8;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

diss_type = 1; % Dissipator type: 0-Poletti, 1-Diehl

%coeff = sqrt(W/J);
coeff = 1;

seed_start = 6;
num_seeds = 1;

init_state_type = 1;
init_state_id = 50;

num_decade_dumps = 50;
begin_decade = -1;
end_decade = 4;
num_decades = end_decade - begin_decade;
num_log_dumps = num_decade_dumps * num_decades + 1;

gamma_id = 1;
data_st = 0;
gammas = 0;

seed = seed_start;

for g = 0.0001:0.0001:0.0009
    local_path = sprintf('results/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.4f/init_type_%d/init_id_%d', diss_type, Nc, W, U, J, g, init_state_type, init_state_id);
    file_name = sprintf('%s/seed_%d/stationary_info_Nc(%d)_dt(%d)_W(%0.2f)_U(%0.2f)_J(%0.2f)_gamma(%0.2f)_it(%d)_is(%d)_seed%d.txt', local_path, seed, Nc, diss_type, W, U, J, g, init_state_type, init_state_id, seed);
	file_name = file_name
    data = importdata(file_name);
    
    size(data)
    
    gammas(gamma_id) = g;
    data_st(gamma_id) = data(5);
    
    gamma_id = gamma_id + 1;
end

for g = 0.001:0.001:0.009
    local_path = sprintf('results/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.4f/init_type_%d/init_id_%d', diss_type, Nc, W, U, J, g, init_state_type, init_state_id);
    file_name = sprintf('%s/seed_%d/stationary_info_Nc(%d)_dt(%d)_W(%0.2f)_U(%0.2f)_J(%0.2f)_gamma(%0.2f)_it(%d)_is(%d)_seed%d.txt', local_path, seed, Nc, diss_type, W, U, J, g, init_state_type, init_state_id, seed);
    data = importdata(file_name);
    
    size(data)
    
    gammas(gamma_id) = g;
    data_st(gamma_id) = data(5);
    
    gamma_id = gamma_id + 1;
end

for g = 0.01:0.01:0.09
    local_path = sprintf('results/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.4f/init_type_%d/init_id_%d', diss_type, Nc, W, U, J, g, init_state_type, init_state_id);
    file_name = sprintf('%s/seed_%d/stationary_info_Nc(%d)_dt(%d)_W(%0.2f)_U(%0.2f)_J(%0.2f)_gamma(%0.2f)_it(%d)_is(%d)_seed%d.txt', local_path, seed, Nc, diss_type, W, U, J, g, init_state_type, init_state_id, seed);
    data = importdata(file_name);
    
    size(data)
    
    gammas(gamma_id) = g;
    data_st(gamma_id) = data(5);
    
    gamma_id = gamma_id + 1;
end

for g = 0.1:0.1:1.00
    local_path = sprintf('results/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.4f/init_type_%d/init_id_%d', diss_type, Nc, W, U, J, g, init_state_type, init_state_id);
    file_name = sprintf('%s/seed_%d/stationary_info_Nc(%d)_dt(%d)_W(%0.2f)_U(%0.2f)_J(%0.2f)_gamma(%0.2f)_it(%d)_is(%d)_seed%d.txt', local_path, seed, Nc, diss_type, W, U, J, g, init_state_type, init_state_id, seed);
    data = importdata(file_name);
    
    size(data)
    
    gammas(gamma_id) = g;
    data_st(gamma_id) = data(5);
    
    gamma_id = gamma_id + 1;
end

hLine = plot(gammas, data_st, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$\gamma$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$S_{st}$', 'Interpreter', 'latex');
set(gca,'xscale','log');
hold all;


