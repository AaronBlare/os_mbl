clear all;

Nc = 8;
W = 10.0;
U = 1.0;
J = 1.0;
g = 0.1;
dt = 1;
alpha = 3.0;
et = 0;
seed = 1;
max_num_seeds = 1000;

Np = Nc/2;
Ns = nchoosek(Nc,Np); 

file_name_suffix = sprintf('_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', Nc, dt, alpha, et, W, U, J, g, max_num_seeds, seed);

file_name = sprintf('hamiltonian%s', file_name_suffix);
data = importdata(file_name);

H_c = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        H_c(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

file_name = sprintf('hamiltonian_ev%s', file_name_suffix);
data = importdata(file_name);

H_ev_c = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        H_ev_c(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

file_name = sprintf('hamiltonian_eg%s', file_name_suffix);
data = importdata(file_name);

H_eg_c = zeros(Ns, 1);
for state_id_1 = 1:Ns
    H_eg_c(state_id_1) = data(state_id_1, 1);
end

file_name = sprintf('trans_rates%s', file_name_suffix);
data = importdata(file_name);

tr_c = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        tr_c(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

file_name = sprintf('rho_in_d%s', file_name_suffix);
data = importdata(file_name);

rho_in_d_c = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        rho_in_d_c(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end