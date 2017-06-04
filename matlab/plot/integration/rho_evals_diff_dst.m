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
W = 2.0;
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

states_num = nchoosek(Nc, Nc/2);

seed_start_begin = seed_start;
seed_start_num = 100;

num_int = 100;
rho_evals_diff_begin = 0;
rho_evals_diff_end = 0.0005;
rho_evals_diff_shift = (rho_evals_diff_end - rho_evals_diff_begin) / num_int;
rho_evals_diff_int = zeros(num_int, 1);
rho_evals_diff_pdf = zeros(num_int, 1);
for int_id = 1:num_int
    rho_evals_diff_int(int_id) = rho_evals_diff_begin + int_id * rho_evals_diff_shift - 0.5 * rho_evals_diff_shift;
end
eps = 1.0e-6;

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
        
        for st_id = 1:states_num - 1
            curr_eval_diff = curr_evals(st_id + 1) - curr_evals(st_id);
			
			if (curr_eval_diff < rho_evals_diff_end)
				int_id = floor((curr_eval_diff - rho_evals_diff_begin) * num_int / (rho_evals_diff_end - rho_evals_diff_begin + eps)) + 1;
				rho_evals_diff_pdf(int_id) = rho_evals_diff_pdf(int_id) + 1;
			end
        end
        
    end 
end

rho_evals_diff_pdf = rho_evals_diff_pdf / (seed_num * seed_start_num * (states_num-1) * rho_evals_diff_shift);

hLine = plot(rho_evals_diff_int, rho_evals_diff_pdf, 'LineWidth', 2);
legend(hLine, sprintf('W=%0.4f \\gamma=%0.4f', W, g));
set(gca, 'FontSize', 30);
xlabel('$\sigma$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF(\sigma)$', 'Interpreter', 'latex');
legend('-DynamicLegend');
legend('Location','southeast')

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

savefig(sprintf('%s/rho_evals_diff_dst_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_evals_diff_dst_%s.pdf', home_figures_path, suffix));


