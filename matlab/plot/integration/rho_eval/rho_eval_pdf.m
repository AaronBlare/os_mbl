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

non_inc_count = 0;

bin_begin = 0.00000001;
num_decades = 10;
bin_end = bin_begin * 10.^num_decades;
num_bin_per_decade = 20;

num_bin = num_bin_per_decade * num_decades;

bin_borders = zeros(num_bin + 1, 1);
bin_centers = zeros(num_bin, 1);
bin_pdf = zeros(num_bin, 1);

for bin_id = 1 : num_bin + 1
    bin_borders(bin_id) = bin_begin * 10.^((bin_id - 1) / num_bin_per_decade);
    if (bin_id <= num_bin)
        bin_centers(bin_id) = bin_begin * 10.^((bin_id - 1 + 0.5) / num_bin_per_decade);
    end
end

bin_diff = diff(bin_borders);

num_evals = size(all_evals, 1);

for eval_id = 1:num_evals
    
    curr_eval = all_evals(eval_id);
    
    if ((curr_eval >= bin_begin) && (curr_eval <= bin_end))
        bin_id = floor((log10(curr_eval) - log10(bin_begin)) * num_bin / (log10(bin_end) - log10(bin_begin) + eps)) + 1;
        bin_pdf(bin_id) = bin_pdf(bin_id) + 1;
    else
        non_inc_count = non_inc_count + 1;
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

non_inc_count = non_inc_count
norm = sum(bin_pdf)




figure
	
figure;
hLine = plot(bin_centers, bin_pdf, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$eval$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
set(gca,'XScale','log');
set(gca,'YScale','log');
ylabel('$PDF$', 'Interpreter', 'latex');

savefig(sprintf('%s/rho_eval_inbin_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_eval_inbin_%s.pdf', home_figures_path, suffix));




for bin_id = 1 : num_bin
    bin_pdf(bin_id) = bin_pdf(bin_id) / (norm * bin_diff(bin_id));
end

norm_check = 0;
for bin_id = 1 : num_bin
    norm_check = norm_check + bin_pdf(bin_id) * bin_diff(bin_id);
end
norm_check = norm_check



figure
	
figure;
hLine = plot(bin_centers, bin_pdf, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$eval$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
set(gca,'XScale','log');
set(gca,'YScale','log');
ylabel('$PDF$', 'Interpreter', 'latex');

savefig(sprintf('%s/rho_eval_pdf_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_eval_pdf_%s.pdf', home_figures_path, suffix));




