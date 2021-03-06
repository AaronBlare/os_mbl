clear all;

data_path = '/data/biophys/yusipov/mbl_zero_super_vector/results';
figures_path = '/home/yusipov/mbl_zero_super_vector/figures';

Nc                  = 8;
dissipator_type     = 1;
alpha               = 0;
energy_type         = 0;
border_conditions   = 1;
W                   = 8;
U                   = 1;
J                   = 1;
g                   = 0.1;
seed_start          = 1;
seed_num            = 100;

num_seed_start = 10;

size_total = seed_num * num_seed_start;

W_start = 1.0;
W_shift = 1.0;
W_num = 10;

Ws = zeros(W_num, 1);
data = zeros(W_num, 1);
W_id = 0;

for W = W_start : W_shift : W_start + W_shift * (W_num-1)
    
    W_id = W_id + 1;
    W = W
    Ws(W_id) = W;
    
    avg_data = 0;
    
    for ss = seed_start : seed_num: seed_start + seed_num * (num_seed_start - 1)
       
        for seed = ss : ss + (seed_num-1)
            
            file_name_tree = sprintf('Nc_%d/dt_%d/alpha_%0.4f/et_%d/bc_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/seed_start_%d', ...
                Nc, ...
                dissipator_type, ...
                alpha, ...
                energy_type, ...
                border_conditions, ...
                W, ...
                U, ...
                J, ...
                g, ...
                ss);
            
            file_name_suffix = sprintf('Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ss(%d)_sn(%d)_seed(%d)', ...
                Nc, ...
                dissipator_type, ...
                alpha, ...
                energy_type, ...
                border_conditions, ...
                W, ...
                U, ...
                J, ...
                g, ...
                ss, ...
                seed_num, ...
                seed);
            
            file_name = sprintf('%s/%s/characteristics_%s.txt.txt', ...
                data_path, file_name_tree, file_name_suffix);
            
            local_data = importdata(file_name);
            
            avg_data = avg_data + abs(local_data(1));
        end
    end
    
    avg_data = avg_data / size_total;
    
    data(W_id) = avg_data;
    
end

file_name_suffix = sprintf('Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_bc(%d)_W(var)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_sn(%d)', ...
    Nc, ...
    dissipator_type, ...
    alpha, ...
    energy_type, ...
    border_conditions, ...
    U, ...
    J, ...
    g, ...
    size_total);

figure;
hLine = plot(Ws, data, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$\W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$<S>$', 'Interpreter', 'latex');

savefig(sprintf('%s/evolution_avg_entropy_from_W_%s.fig', ...
    figures_path, file_name_suffix));