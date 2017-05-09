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
num_seeds = 10000;

num_intervals_imb = 100;

imbalances = zeros(num_seeds, 1);

for seed = seed_start : seed_start + (num_seeds - 1)

	seed = seed
	
    path = sprintf('../results/Nc_%d/dt_%d/alpha_%0.4f/et_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/max_num_seeds_%d/seed_%d', ...
        Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
    
    file_name = sprintf('%s/characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', ...
        path, Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);
		
    data = importdata(file_name);
    imbalances(seed - seed_start + 1) = data(2);
end

min_imb = min(imbalances);
max_imb = max(imbalances);

shift_imb = (max_imb - min_imb) / num_intervals_imb;
intervals_imb = zeros(num_intervals_imb, 1);
pdf_imb = zeros(num_intervals_imb, 1);

for i = 1:num_intervals_imb
    intervals_imb(i) = min_imb + i * shift_imb - 0.5 * shift_imb;
end

for i = 1:num_seeds
    
    current_imb = imbalances(i);
    
    if ((current_imb >= min_imb) && (current_imb <= max_imb))
        id_imb = floor((current_imb - min_imb) * num_intervals_imb / (max_imb - min_imb + 0.000001)) + 1;
        pdf_imb(id_imb) = pdf_imb(id_imb) + 1;
    end
end

for i = 1:num_intervals_imb
    pdf_imb(i) = pdf_imb(i) / (num_seeds * shift_imb);
end


figure;
hLine = plot(intervals_imb, pdf_imb, '-o', 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$I$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF(I)$', 'Interpreter', 'latex');

savefig(sprintf('pdf_imbalance_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_ss(%d)_sn(%d).fig', ...
	Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed_start, num_seeds));