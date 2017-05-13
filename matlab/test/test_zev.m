clear all;

filename = '../source/zev/config.txt';
input_data = importdata(filename);

Nc                  = input_data(1);
dissipator_type     = input_data(2);
alpha               = input_data(3);		
energy_type         = input_data(4);
border_conditions   = input_data(5);
W                   = input_data(6);
U                   = input_data(7);
J                   = input_data(8);
g                   = input_data(9);
seed_start          = input_data(10);
seed_num            = input_data(11);
save_type           = input_data(12);
save_super_operator = input_data(13);
dump_type           = input_data(14);

max_num_seeds = 1000;
seed_cpp = 1;

Np = Nc/2;
Ns = nchoosek(Nc,Np); 

path_matlab = '../../data/zev/matlab/';
path_cpp = '../../data/zev/cpp/';



suffix_matlab = sprintf('_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ss(%d)_sn(%d)_seed(%d).txt', ...
    Nc, ...
    dissipator_type, ...
    alpha, ...
    energy_type, ...
    border_conditions, ...
    W, ...
    U, ...
    J, ...
    g, ...
    seed_start, ...
    seed_num, ...
    seed_start);

suffix_cpp = sprintf('_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', ...
    Nc, ...
    dissipator_type, ...
    alpha, ...
    energy_type, ...
    W, ...
    U, ...
    J, ...
    g, ...
    max_num_seeds, ...
    seed_cpp);

file_name = sprintf('%shamiltonian%s', path_matlab, suffix_matlab);
data = importdata(file_name);
mtx_matlab = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

data = 0;

file_name = sprintf('%shamiltonian%s', path_cpp, suffix_cpp);
data = importdata(file_name);
mtx_cpp = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1);
    end
end

H_diff = max(max(abs(mtx_matlab - mtx_cpp)))





% for dissipator_id = 1:Nc-1
%     
%     file_name = sprintf('%sdissipator_%d%s', path_matlab, dissipator_id, suffix_matlab);
%     data = importdata(file_name);
%     mtx_matlab = zeros(Ns, Ns);
%     for state_id_1 = 1:Ns
%         for state_id_2 = 1:Ns
%             mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
%         end
%     end
%     
%     data = 0;
%     
%     file_name = sprintf('%sdissipator_%d%s', path_cpp, dissipator_id, suffix_cpp);
%     data = importdata(file_name);
%     mtx_cpp = zeros(Ns, Ns);
%     for state_id_1 = 1:Ns
%         for state_id_2 = 1:Ns
%             mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
%         end
%     end
%     
%     dissipator_diff = max(max(abs(mtx_matlab - mtx_cpp)))
%     
%     file_name = sprintf('%sdissipator_conj_%d%s', path_matlab, dissipator_id, suffix_matlab);
%     data = importdata(file_name);
%     mtx_matlab = zeros(Ns, Ns);
%     for state_id_1 = 1:Ns
%         for state_id_2 = 1:Ns
%             mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
%         end
%     end
%     
%     data = 0;
%     
%     file_name = sprintf('%sdissipator_conj_%d%s', path_cpp, dissipator_id, suffix_cpp);
%     data = importdata(file_name);
%     mtx_cpp = zeros(Ns, Ns);
%     for state_id_1 = 1:Ns
%         for state_id_2 = 1:Ns
%             mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
%         end
%     end
%     
%     dissipator_conj_diff = max(max(abs(mtx_matlab - mtx_cpp)))
%     
%     file_name = sprintf('%sdissipator_mult_%d%s', path_matlab, dissipator_id, suffix_matlab);
%     data = importdata(file_name);
%     mtx_matlab = zeros(Ns, Ns);
%     for state_id_1 = 1:Ns
%         for state_id_2 = 1:Ns
%             mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
%         end
%     end
%     
%     data = 0;
%     
%     file_name = sprintf('%sdissipator_mult_%d%s', path_cpp, dissipator_id, suffix_cpp);
%     data = importdata(file_name);
%     mtx_cpp = zeros(Ns, Ns);
%     for state_id_1 = 1:Ns
%         for state_id_2 = 1:Ns
%             mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
%         end
%     end
%     
%     dissipator_mult_diff = max(max(abs(mtx_matlab - mtx_cpp)))
%     
% end










file_name = sprintf('%slindbladian%s', path_matlab, suffix_matlab);
data = importdata(file_name);
mtx_matlab = zeros(Ns, Ns);
for state_id_1 = 1:Ns*Ns
    for state_id_2 = 1:Ns*Ns
        mtx_matlab(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
    end
end

data = 0;

file_name = sprintf('%slindbladian%s', path_cpp, suffix_cpp);
data = importdata(file_name);
mtx_cpp = zeros(Ns, Ns);
for state_id_1 = 1:Ns*Ns
    for state_id_2 = 1:Ns*Ns
        mtx_cpp(state_id_1, state_id_2) = data((state_id_1-1) * Ns + state_id_2, 1) + sqrt(-1) * data((state_id_1-1) * Ns + state_id_2, 2);
    end
end

lindbladian_diff = max(max(abs(mtx_matlab - mtx_cpp)))