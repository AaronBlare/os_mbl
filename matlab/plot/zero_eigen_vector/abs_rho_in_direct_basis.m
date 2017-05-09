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

seed = 6;

path = sprintf('../results/dt_%d/alpha_%0.4f/et_%d/Ns_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/seed_%d', ...
    diss_type,alpha, energy_type, Nc, W, U, J, g, seed);

file_name = sprintf('%s/rho_in_direct_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', ...
    path, Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
	
file_name = file_name

data = importdata(file_name);

abs_rho = zeros(Ns, Ns);
states = zeros(Ns, 1);
for state_id_1 = 1:Ns
    
    states(state_id_1) = state_id_1;
    
    for state_id_2 = 1:Ns
        abs_rho(state_id_1, state_id_2) = abs( data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2));
    end
end

figure;
hLine = imagesc(states, states, abs_rho);
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$m$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '|\rho_{n,m}|');
set(gca,'YDir','normal');

savefig(sprintf('abs_rho_in_direct_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).fig', ...
	Nc, diss_type, alpha, energy_type, W, U, J, g, seed));