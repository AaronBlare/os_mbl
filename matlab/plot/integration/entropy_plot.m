clear all;

W = 12; % disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

Nc = 8;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

diss_type = 1; % Dissipator type: 0-Poletti, 1-Diehl

%coeff = sqrt(W/J);
coeff = 1;

seed_start = 2;
num_seeds = 1;

init_state_type = 1;
init_state_id = 50;

num_decade_dumps = 50;
begin_decade = -1;
end_decade = 4;
num_decades = end_decade - begin_decade;
num_log_dumps = num_decade_dumps * num_decades + 1;

energy_type = 1;

local_path = sprintf('results/dt_%d/et_%d/Ns_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/init_type_%d/init_id_%d', diss_type, energy_type, Nc, W, U, J, g, init_state_type, init_state_id);

file_name = sprintf('%s/seed_%d/times_Nc(%d)_dt(%d)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_it(%d)_is(%d).txt', local_path, seed_start, Nc, diss_type, energy_type, W, U, J, g, init_state_type, init_state_id);
file_name = file_name
log_times = importdata(file_name);
log_times = log_times(2:end);

size(log_times, 1)

entropies = zeros(num_log_dumps, 1);

for seed = seed_start : seed_start + (num_seeds - 1)

	seed = seed
     
    file_name = sprintf('%s/seed_%d/entropy_Nc(%d)_dt(%d)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_it(%d)_is(%d)_seed%d.txt', local_path, seed, Nc, diss_type, energy_type, W, U, J, g, init_state_type, init_state_id, seed);
    local_entropies = importdata(file_name);
	local_entropies = local_entropies(1+1 : num_log_dumps+1);
    
    for dump_id = 1:num_log_dumps
        entropies(dump_id) = entropies(dump_id) + local_entropies(dump_id) / num_seeds;
    end
    
end

hLine = plot(log_times * g, coeff * entropies, 'LineWidth', 2);
%hLine = plot(log_times, coeff * entropies, 'LineWidth', 2);
legend(hLine, sprintf('W=%0.4f \\gamma=%0.4f', W, g));
set(gca, 'FontSize', 30);
xlabel('$\gamma t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$S(\gamma t)$', 'Interpreter', 'latex');
set(gca,'xscale','log');

xlim([10^begin_decade 10^end_decade])
legend('-DynamicLegend');
legend('Location','southeast')

hold all;


