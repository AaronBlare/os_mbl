clear all;

home_figures_path = '/home/yusipov/Work/os_mbl/figures/int/matlab';

data_path = '/data/biophys/yusipov/os_mbl/int/matlab/';
prefix = 'characteristics';
data_path = sprintf('%s%s', data_path, prefix);

Nc = 8;
diss_type = 1;
diss_phase = 0.0;
energy_type = 0;
periodic_bc = 0;
W = 10.0;
U = 1.0;
J = 1.0;
g = 0.1;
seed_start = 1;
seed_num = 100;
is_int = 0;
int_ist = 0;
int_isi = 50;
int_dt = 1;
int_db = 0.1;
int_de = 10000;
int_dn = 250;
is_zev = 1;
is_zev_check = 1;
is_eg_stat = 1;
is_save_vec = 1;
is_save_mtx = 0;
fs_type = 0;

seed_start_begin = seed_start;
seed_start_num = 100;

states_num = nchoosek(Nc, Nc/2);
    
all_evals = [];

for seed_start = seed_start_begin : seed_num : seed_start_begin + seed_num * (seed_start_num - 1)
    
    for seed = seed_start : seed_start + (seed_num-1)
	
		seed = seed
        
        curr_path = sprintf('%s/Nc_%d/dt_%d/dp_%0.4f/et_%d/bc_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/seed_start%d', ...
            data_path, ...
            Nc, ...
            diss_type, ...
            diss_phase, ...
            energy_type, ...
            periodic_bc, ...
            W, ...
            U, ...
            J, ...
            g, ...
            seed_start);
        
        suffix = sprintf('zev_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
            Nc, ...
            diss_type, ...
            diss_phase, ...
            energy_type, ...
            periodic_bc, ...
            W, ...
            U, ...
            J, ...
            g, ...
            seed);
        
        file_name = sprintf('%s/rho_evals_%s', curr_path, suffix);
        
        curr_evals = importdata(file_name);
        curr_evals = sort(curr_evals);
        
        all_evals = vertcat(all_evals, curr_evals);
        
    end
end

all_evals = sort(all_evals);

num_int = 1000;
evals_begin = min(all_evals)
evals_end = max(all_evals)
evals_shift = (evals_end - evals_begin) / num_int;
evals_int = zeros(num_int, 1);
evals_hist = zeros(num_int, 1);
for int_id = 1:num_int
    evals_int(int_id) = evals_begin + int_id * evals_shift - 0.5 * evals_shift;
end
eps = 1.0e-8;

num_evals = size(all_evals, 1);

for eval_id = 1:num_evals
    
    curr_eval = all_evals(eval_id);
    
    if (curr_eval <= evals_end && curr_eval >= evals_begin)
        int_id = floor((curr_eval - evals_begin) * num_int / (evals_end - evals_begin + eps)) + 1;
        evals_hist(int_id) = evals_hist(int_id) + 1;
    end
end

suffix = sprintf('Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(var)', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    periodic_bc, ...
	W, ...
    U, ...
    J, ...
    g);

figure
	
hLine = plot(evals_int, evals_hist, 'LineWidth', 2);
legend(hLine, sprintf('W=%0.4f', W));
set(gca, 'FontSize', 30);
xlabel('$eval$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
legend('-DynamicLegend');
legend('Location','southeast')

savefig(sprintf('%s/rho_eval_hist_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_eval_diff_pdf_%s.pdf', home_figures_path, suffix));




