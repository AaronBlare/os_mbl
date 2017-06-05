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

W_begin = 0.1;
W_shift = 0.1;
W_num = 100;

seed_start_begin = seed_start;
seed_start_num = 100;

states_num = nchoosek(Nc, Nc/2);

num_int = 500;
evals_diff_begin = 0;
evals_diff_end = 0.00025;
evals_diff_shift = (evals_diff_end - evals_diff_begin) / num_int;
evals_diff_int = zeros(num_int, 1);
evals_diff_pdf = zeros(num_int, 1);

imb_evo = zeros(W_num, num_int);

Ws = zeros(W_num, 1);
for W_id = 1:W_num
    
    W = W_begin + (W_id-1) * W_shift
    Ws(W_id) = W;
    
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
	
	file_name = file_name
    
    imb_data = importdata(file_name);
    
    imb_int = imb_data(:,1);
    imb_evo(W_id, :) = imb_data(:,2);
    
end

suffix = sprintf('Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(var)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(var)', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
    periodic_bc, ...
    U, ...
    J, ...
    g);

hLine = imagesc(Ws, imb_int, imb_evo');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\sigma$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF(\sigma)', 'FontSize', 33);
set(gca,'YDir','normal');
hold all;

savefig(sprintf('%s/rho_eval_diff_pdf_evo_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_eval_diff_pdf_evo_%s.pdf', home_figures_path, suffix));


hLine = imagesc(Ws, imb_int, log10(imb_evo' + 1.0e-6));
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$I$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'log_{10}PDF(I)', 'FontSize', 33);
set(gca,'YDir','normal');
hold all;

savefig(sprintf('%s/rho_eval_diff_pdf_evo_log_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(gcf, 'renderer','painters');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_eval_diff_pdf_evo_log_%s.pdf', home_figures_path, suffix));


