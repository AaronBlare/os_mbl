clear all;

Nc = 14;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

dt = 1; 		% Dissipator type: 0-Poletti, 1-Diehl
alpha = 0; 			% Dissipator Diehl phase 
et = 0;   	% 0 if regular, 1 if zero mean

W = 8; 	% disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

max_num_seeds = 1000000;
seed_start = 1;
num_seeds = 1000;

W_start = 1;
num_Ws = 10;

Ws = zeros(num_Ws, 1);
abs_ent = zeros(num_Ws, 1);

for W = W_start : W_start + (num_Ws - 1)
    
    W = W
    Ws(W - W_start + 1) = W;
    
    local_ent = 0.0;
    for seed = seed_start : seed_start + (num_seeds - 1)
	
		seed = seed
        
        path = sprintf('../results/Nc_%d/dt_%d/alpha_%0.4f/et_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/max_num_seeds_%d/seed_%d', ...
            Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
        
        file_name = sprintf('%s/characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', ...
            path, Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
        
        data = importdata(file_name);
        local_ent = local_ent + data(1);
    end
    
    local_ent = local_ent / num_seeds;
    
    abs_ent(W - W_start + 1) = local_ent;
    
end

figure;
hLine = plot(Ws, abs_ent, '-o', 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$<S>$', 'Interpreter', 'latex');

savefig(sprintf('evo_entropy_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(var)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_ss(%d)_sn(%d).fig', ...
	Nc, dt, alpha, et, U, J, g, max_num_seeds, seed_start, num_seeds));