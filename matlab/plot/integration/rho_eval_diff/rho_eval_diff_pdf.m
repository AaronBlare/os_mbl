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
W = 10;
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

W = 10.0;

seed_start_begin = seed_start;
seed_start_num = 100;

states_num = nchoosek(Nc, Nc/2);

num_int = 500;
evals_diff_begin = 0;
evals_diff_end = 0.00025;
evals_diff_shift = (evals_diff_end - evals_diff_begin) / num_int;
evals_diff_int = zeros(num_int, 1);
evals_diff_pdf = zeros(num_int, 1);

curr_path = sprintf('%s/Nc_%d/dt_%d/dp_%0.4f/et_%d/bc_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/rho_eval_diff_pdf', ...
    data_path, ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    periodic_bc, ...
    W, ...
    U, ...
    J, ...
    g);

suffix = sprintf('_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(var).txt', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    periodic_bc, ...
    W, ...
    U, ...
    J, ...
    g);

file_name = sprintf('%s/rho_evals_diff_pdf%s', curr_path, suffix);

imb_data = importdata(file_name);

imb_int = imb_data(:,1);
imb_evo = imb_data(:,2);

norm = sum(imb_evo) * 1.0 / 5000.0

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
	
hLine = plot(imb_int, imb_evo, 'LineWidth', 2);
legend(hLine, sprintf('W=%0.4f', W));
set(gca, 'FontSize', 30);
xlabel('$\sigma$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF(\sigma)$', 'Interpreter', 'latex');
legend('-DynamicLegend');
legend('Location','southeast')

savefig(sprintf('%s/rho_eval_diff_pdf_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_eval_diff_pdf_%s.pdf', home_figures_path, suffix));

file_name = sprintf('rho_eval_diff_pdf_%s.txt', suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:size(imb_int, 1)
	fprintf(file_id, '%0.16e %0.16e\n', imb_int(dump_id), imb_evo(dump_id));
end
fclose(file_id);
