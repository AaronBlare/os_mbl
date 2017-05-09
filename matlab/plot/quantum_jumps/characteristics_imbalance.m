clear all;

imag1 = sqrt(-1);


W = 10; % disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

seed_start = 6;
seed_num = 1;

dt = 1; % dissipator type: 1-Poletti, 2-Diehl

Nc = 8;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states
num = precalc_states (Nc, Np);

period = 0.1;
num_periods = 100000;
init_state_id = 49;
num_dumps = 251;
dump_type = 1;
num_periods_in_trans_proc = 0;
num_trajectories = 10000;
rnd_cur = 0;

imbalance = zeros(num_dumps, 1);
time_dumps = zeros(num_dumps, 1);

for seed = seed_start : seed_start + (seed_num-1)

	seed = seed
    
    if seed == seed_start
        file_name = sprintf('../../results/period_%0.2f/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.2f/is_%d/seed_%d/dump_%d/rnd_%d/periods.txt', period, dt, Ns, W, U, J, g, init_state_id, seed, dump_type, rnd_cur);
        periods_data = importdata(file_name);
        time_dumps = periods_data(1 : num_dumps);
		periods = time_dumps;
        
        time_dumps = time_dumps * period;
    end
    
    file_name = sprintf('../../results/period_%0.2f/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.2f/is_%d/seed_%d/dump_%d/rnd_%d/characteristics.txt', period, dt, Ns, W, U, J, g, init_state_id, seed, dump_type, rnd_cur);
	characteristics_data = importdata(file_name);
    
    current_imbalance = characteristics_data(:, 2);
	
	size(current_imbalance)
    
    imbalance = imbalance + current_imbalance / seed_num;
    
end

hLine = plot(time_dumps * g, -log(imbalance));
legend(hLine, sprintf('W=%0.2f \\gamma=%0.2f', W, g));
set(gca, 'FontSize', 30);
xlabel('$\gamma t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$-\log I(\gamma t)$', 'Interpreter', 'latex');
set(gca,'xscale','log');
set(gca,'yscale','log');
hold all;