clear all;

imag1 = sqrt(-1);


W = 10; % disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

seed_start = 1;
seed_num = 1;

dt = 1; % dissipator type: 1-Poletti, 2-Diehl

Nc = 8;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states
num = precalc_states (Nc, Np);

period = 0.1;
num_periods = 100000;
init_state_id = 49;
num_dumps = 51;
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
    
    for i = 1:num_dumps
	
		file_name = sprintf('../../results/period_%0.2f/dt_%d/Ns_%d/W_%0.2f/U_%0.2f/J_%0.2f/g_%0.2f/is_%d/seed_%d/dump_%d/rnd_%d/rho_period_%d.txt', period, dt, Ns, W, U, J, g, init_state_id, seed, dump_type, rnd_cur, periods(i));
		rho_data = importdata(file_name);
	
		rho = zeros(Ns);
        for state_1 = 1:Ns
            for state_2 = 1:Ns
                rho(state_1, state_2) = rho_data((state_1-1) * Ns + state_2, 1) + imag1 * rho_data((state_1-1) * Ns + state_2, 2);
            end
        end
        
        norm = sum(abs(diag(rho)));
		
		i = i
		norm_diff = abs(1.0 - norm)
        
        n_part=zeros(1,Nc);
        for k = 1:Ns
            n_part = n_part + real(rho(k,k)) * ( dec2bin(idtox(k),Nc) == '1' );
        end
        
        curr_imbalance = sum(n_part(1:2:Nc-1) - n_part(2:2:Nc)) / sum(n_part);
        
        imbalance(i) = imbalance(i) + curr_imbalance / seed_num;   
    end    
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

%file_name = sprintf('entropy_Nc(%d)_d(%d)_W(%0.1f)_U(%0.1f)_J(%0.1f)_gamma(%0.1f).fig', Nc, dt, W, U, J, g);
%savefig(file_name)