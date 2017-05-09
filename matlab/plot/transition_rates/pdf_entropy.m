clear all;

Nc = 14;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

dt = 1; 		% Dissipator type: 0-Poletti, 1-Diehl
alpha = 0; 			% Dissipator Diehl phase 
et = 0;   	% 0 if regular, 1 if zero mean

W = 4; 	% disorder
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

max_num_seeds = 1000000;
seed_start = 1;
num_seeds = 10000;

num_intervals_ent = 100;

entropies = zeros(num_seeds, 1);

for seed = seed_start : seed_start + (num_seeds - 1)

	seed = seed
	
    path = sprintf('../results/Nc_%d/dt_%d/alpha_%0.4f/et_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/max_num_seeds_%d/seed_%d', ...
        Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
    
    file_name = sprintf('%s/characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', ...
        path, Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
		
    data = importdata(file_name);
    entropies(seed - seed_start + 1) = data(1);
end

min_ent = min(entropies);
max_ent = max(entropies);

shift_ent = (max_ent - min_ent) / num_intervals_ent;
intervals_ent = zeros(num_intervals_ent, 1);
pdf_ent = zeros(num_intervals_ent, 1);

for i = 1:num_intervals_ent
    intervals_ent(i) = min_ent + i * shift_ent - 0.5 * shift_ent;
end

for i = 1:num_seeds
    
    current_ent = entropies(i);
    
    if ((current_ent >= min_ent) && (current_ent <= max_ent))
        id_ent = floor((current_ent - min_ent) * num_intervals_ent / (max_ent - min_ent + 0.000001)) + 1;
        pdf_ent(id_ent) = pdf_ent(id_ent) + 1;
    end
end

for i = 1:num_intervals_ent
    pdf_ent(i) = pdf_ent(i) / (num_seeds * shift_ent);
end


figure;
hLine = plot(intervals_ent, pdf_ent, '-o', 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$S$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF(S)$', 'Interpreter', 'latex');

savefig(sprintf('pdf_entropy_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_ss(%d)_sn(%d).fig', ...
	Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed_start, num_seeds));