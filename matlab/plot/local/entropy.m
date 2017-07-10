clear all;

task = 2;
Nc = 10;
dt = 1;
dp = 0.0;
et = 0;
bc = 0; 
W = 1.0;
U = 1.0;
J = 1.0;
g = 0.1;

data_path = '../../../data/cluster/unn';

Np = Nc/2;
Ns = nchoosek(Nc, Np);

num_seeds = 100;

ntd = 102;

ee = zeros(num_seeds, ntd);

for seed = 1:num_seeds
    
    fn_path = sprintf('task_%d/Nc_%d/dt_%d/dp_%0.4f/et_%d/bc_%d/W_%0.4f/U_%0.4f/J_%0.4f/g_%0.4f/seed_%d', ...
        task, ...
        Nc, ...
        dt, ...
        dp, ...
        et, ...
        bc, ...
        W, ...
        U, ...
        J, ...
        g, ...
        seed);
    
    fn_suffix = sprintf('f_int_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
        Nc, ...
        dt, ...
        dp, ...
        et, ...
        bc, ...
        W, ...
        U, ...
        J, ...
        g, ...
        seed);
    
    if seed == 1
        fn = sprintf('%s/%s/times_%s', data_path, fn_path, fn_suffix);
        times = importdata(fn);
    end
    
   fn = sprintf('%s/%s/entropy_%s', data_path, fn_path, fn_suffix);
    curr_ee = importdata(fn);
    
    ee(seed, :) = curr_ee(:);
end

ee_mean = zeros(ntd, 1);
ee_m2 = zeros(ntd, 1);

for dump_id = 1:ntd
    for seed = 1:num_seeds
        ee_mean(dump_id) = ee_mean(dump_id) + ee(seed, dump_id);
    end
    ee_mean(dump_id) = ee_mean(dump_id) / num_seeds;
end

for dump_id = 1:ntd
    for seed = 1:num_seeds
        ee_m2(dump_id) = (ee(seed, dump_id) - ee_mean(dump_id)).^2;
    end
    ee_m2(dump_id) = ee_m2(dump_id) / num_seeds;
end

ee_std = sqrt(ee_m2);

hLine = plot(times, ee_mean, 'LineWidth', 2);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
set(gca, 'XScale', 'log');
hold all;

hLine = plot(times, ee_mean + ee_std);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
set(gca, 'XScale', 'log');
hold all;

hLine = plot(times, ee_mean - ee_std);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
set(gca, 'XScale', 'log');
hold all;
