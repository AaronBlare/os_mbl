clear all;

W_start = 1;
W_shift = 1;
W_num = 15;
U = 1;  % interaction
J = 1;  % hopping
g = 0.1;  % gamma;

Nc = 8;                 % number of cells
Np = Nc/2;              % number of particles
Ns = nchoosek(Nc, Np);  % number of states

diss_type = 1; % Dissipator type: 0-Poletti, 1-Diehl

%coeff = sqrt(W/J);
coeff = 1;

seed_start = 10;
num_seeds = 1;

init_state_type = 1;
init_state_id = 50;

num_decade_dumps = 50;
begin_decade = -1;
end_decade = 4;
num_decades = end_decade - begin_decade;
num_log_dumps = num_decade_dumps * num_decades + 1;

energy_type = 0;

Ws = zeros(W_num, 1);
st_entropies_avg = zeros(W_num, 1);
for W_id = 1 : W_num
    W = W_start + W_shift * (W_id - 1);
	W=W
    Ws(W_id) = W;
    
    local_path = sprintf('results/dt_%d/et_%d/Ns_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/init_type_%d/init_id_%d', diss_type, energy_type, Nc, W, U, J, g, init_state_type, init_state_id);
    
    st_entropy_avg = 0;
    
    for seed = seed_start : seed_start + (num_seeds - 1)
        
        file_name = sprintf('%s/seed_%d/stationary_info_Nc(%d)_dt(%d)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_it(%d)_is(%d)_seed%d.txt', local_path, seed, Nc, diss_type, energy_type, W, U, J, g, init_state_type, init_state_id, seed);
        data = importdata(file_name);
        
        st_entropy_avg  = st_entropy_avg + data(5) / num_seeds;
    end
    
    st_entropies_avg(W_id) = st_entropy_avg;
end

Ws

st_entropies_avg

hLine = plot(Ws, st_entropies_avg, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$entropy_{st}$', 'Interpreter', 'latex');
hold all;

