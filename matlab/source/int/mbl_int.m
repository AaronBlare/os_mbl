clear all;

tic

filename = 'config.txt';
input_data = importdata(filename);

Nc                  = input_data(1);
diss_type           = input_data(2);
diss_phase          = input_data(3);
energy_type         = input_data(4);
periodic_bc 		= input_data(5);
W                   = input_data(6);
U                   = input_data(7);
J                   = input_data(8);
g                   = input_data(9);
seed                = input_data(10);
init_state_type     = input_data(11);
init_state_id       = input_data(12);
dump_type           = input_data(13);
begin_dump          = input_data(14);
end_dump            = input_data(15);
num_dumps           = input_data(16);
save_type           = input_data(17);
file_system_type    = input_data(18);

if file_system_type == 1
    data_path = '../../../data/int/matlab/';
elseif file_system_type == 0
    data_path = '';
else
   error('Error: wrong dump_type');
end

Np = Nc/2;
Ns = nchoosek(Nc, Np);
num = precalc_states (Nc, Np, periodic_bc);

file_name_suffix = sprintf('Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ist(%d)_iss(%d)_dump_type(%d)_seed(%d).txt', ...
    Nc, ...
    diss_type, ...
    diss_phase, ...
    energy_type, ...
	periodic_bc, ...
    W, ...
    U, ...
    J, ...
    g, ...
    init_state_type, ...
    init_state_id, ...
    dump_type, ...
    seed);

dump_times = zeros(num_dumps + 1, 1);
if (dump_type == 0)
    dump_shift = (end_dump - begin_dump) / num_dumps;
    for dump_id = 1:num_dumps+1
        dump_times(dump_id) = begin_dump + (dump_id - 1) * dump_shift;
    end
elseif (dump_type == 1)
    begin_decade = log10(begin_dump);
    end_decade = log10(end_dump);
    num_decades = end_decade - begin_decade;
    num_decade_dumps = num_dumps / num_decades;
    for dump_id = 1:num_dumps+1
        dump_times(dump_id) = power(10, begin_decade) * power(10, (1.0/num_decade_dumps) * (dump_id - 1));
    end
else
    error('Error: wrong dump_type');
end

dump_times = vertcat(0, dump_times);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% disorder
%%%%%%%%%%%%%%
Hd = zeros(Ns);
rng(seed)
E = 2 * rand(Nc,1) - 1;

if (energy_type == 1)
    E = E - sum(E)/Nc;
end

for s_id = 1:Ns
    Hd(s_id,s_id) = ( dec2bin(idtox(s_id), Nc) == '1' ) * E;
end

file_name = sprintf('%srandom_energies_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for c_id = 1:Nc
    fprintf(file_id, '%0.18e\n', E(c_id));
end
fclose(file_id);

%%%%%%%%%%%%%%%%%%%%%
% Interaction
%%%%%%%%%%%%%%%%%%%%%

Hi = zeros(Ns);
for s_id_1 = 1:Ns
    
    shifted = bitshift(idtox(s_id_1), 1);
    
    if (periodic_bc == 1)
        if (idtox(s_id_1) >= 2^(Nc-1))
            shifted = shifted + 1;
        end
    end
    
    Hi(s_id_1, s_id_1) = sum( dec2bin(bitand(idtox(s_id_1), shifted)) == '1' );
    
end

%%%%%%%%%%%%%%%%%%%%%
% Hopping
%%%%%%%%%%%%%%%%%%%%%
Hh = zeros(Ns);
for s_id_1 = 1:Ns
    for s_id_2 = 1:Ns
        Hh(s_id_1,s_id_2) = is_adjacent(idtox(s_id_1), idtox(s_id_2));
    end
end

%%%%%%%%%%%%%%%%%%%%%
% Total
%%%%%%%%%%%%%%%%%%%%%
H= -J*Hh + U*Hi + 2.0*W*Hd;

if(save_type >= 1)
    file_name = sprintf('%shamiltonian_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for s_id_1 = 1:Ns
        for s_id_2 = 1:Ns
            fprintf(file_id, '%0.18e\n', H(s_id_1, s_id_2));
        end
    end
    fclose(file_id);
end

[Ev,Eg] = eig(H); % Anderson modes

if(save_type >= 1)
    file_name = sprintf('%seg_hamiltonian_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for s_id = 1:Ns
        fprintf(file_id, '%0.18e\n', Eg(s_id, s_id));
    end
    fclose(file_id);
end

if(save_type >= 2)
    file_name = sprintf('%sev_hamiltonian_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for s_id_1 = 1:Ns
        for s_id_2 = 1:Ns
            fprintf(file_id, '%0.18e\n', Ev(s_id_1, s_id_2));
        end
    end
    fclose(file_id);
end

%%%%%%%%%%%%%%%%%%%
% Creating matrix of right part P*Rho=0
%%%%%%%%%%%%%%%%%%%
lndbldn = zeros((Ns)^2);
lndbldn = -sqrt(-1) * (kron(eye(Ns),H) - kron(transpose(H),eye(Ns)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (diss_type == 1)
    
    for diss_id = 1:Nc-1
        diss = zeros(Ns);
        for s_id_1 = 1:Ns
            diss(s_id_1,s_id_1) = bitget(idtox(s_id_1),diss_id) - bitget(idtox(s_id_1),diss_id+1);
            for s_id_2 = 1:Ns
                if(is_adjacent(idtox(s_id_1), idtox(s_id_2)))
                    hop = 1 + Nc - find(dec2bin(bitxor(idtox(s_id_1),idtox(s_id_2)),Nc)=='1');
                    if(hop(2) == diss_id)
                        if((bitget(idtox(s_id_1),diss_id)))
                            diss(s_id_2,s_id_1) = exp(sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = -exp(-sqrt(-1)*diss_phase);
                        else
                            diss(s_id_2,s_id_1) = -exp(-sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = exp(sqrt(-1)*diss_phase);
                        end
                    end
                end
            end
        end
        
        lndbldn = lndbldn + ...
            g * 0.5 *(2 * kron(eye(Ns),diss) * kron(transpose(diss'),eye(Ns)) - ...
            kron(transpose(diss'*diss),eye(Ns)) - kron(eye(Ns),diss'*diss));
        
    end
    
    if (periodic_bc == 1)
        
        diss = zeros(Ns);
        for s_id_1 = 1:Ns
            diss(s_id_1, s_id_1) = bitget(idtox(s_id_1), Nc) - bitget(idtox(s_id_1), 1);
            for s_id_2 = 1:Ns
                if(is_adjacent(idtox(s_id_1), idtox(s_id_2)))
                    hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(s_id_1), idtox(s_id_2)),Nc) == '1' );
                    if(hopping_ids(2) == Nc)
                        if(bitget(idtox(s_id_1), Nc))
                            diss(s_id_2,s_id_1) = exp(sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = -exp(-sqrt(-1)*diss_phase);
                        else
                            diss(s_id_2,s_id_1) = -exp(-sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = exp(sqrt(-1)*diss_phase);
                        end
                    end
                end
            end
        end
        
       lndbldn = lndbldn + ...
            g * 0.5 *(2 * kron(eye(Ns),diss) * kron(transpose(diss'),eye(Ns)) - ...
            kron(transpose(diss'*diss),eye(Ns)) - kron(eye(Ns),diss'*diss));
    end
    
elseif(diss_type == 0)
    
    for diss_id = 1:Nc
        diss = zeros(Ns);
        for s_id=1:Ns
            diss(s_id,s_id) = bitget(idtox(s_id),diss_id);
        end
        
        lndbldn = lndbldn + ...
            g * 0.5 *(2 * kron(eye(Ns),diss) * kron(transpose(diss'),eye(Ns)) - ...
            kron(transpose(diss'*diss),eye(Ns)) - kron(eye(Ns),diss'*diss));
    end
else
    error('Error: wrong dissipator_type');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

lndbldn_sprs = sparse(lndbldn);
lndbldn = 0;

start_rho = zeros(Ns * Ns, 1);
if init_state_type == 0
    start_rho((init_state_id - 1) * (Ns) + init_state_id) = 1.0;
elseif init_state_type == 1
    for s_id = 1:Ns
        start_rho((s_id - 1) * (Ns) + s_id) = 1.0 / Ns;
    end
else
    error('Error: wrong init_state_type');
end

tic
[times, rho_dumps] = ode45(@(times, rho_dumps) right_part(times, rho_dumps, lndbldn_sprs), dump_times, start_rho);
toc

total_num_dumps = size(times, 1);

file_name = sprintf('%stimes_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', times(dump_id));
end
fclose(file_id);

local_entropies         = zeros(total_num_dumps, 1);
local_entropies_and     = zeros(total_num_dumps, 1);
local_imbalances        = zeros(total_num_dumps, 1);
local_imbalances_and    = zeros(total_num_dumps, 1);
local_iprs              = zeros(total_num_dumps, 1);
local_iprs_and          = zeros(total_num_dumps, 1);

final_rho = zeros(Ns, Ns);
final_rho_and = zeros(Ns, Ns);

for dump_id = 1:total_num_dumps
    
    current_rho_vec = rho_dumps(dump_id, :);
    
    current_rho_mat = zeros(Ns,Ns);
    for s_id_1 = 1:Ns
        for s_id_2 = 1:Ns
            current_rho_mat(s_id_1, s_id_2) = current_rho_vec((s_id_1-1)*(Ns) + s_id_2);
        end
    end
    
    if (dump_id == total_num_dumps)
        final_rho = current_rho_mat;
    end
    
    current_rho_mat_and = Ev' * current_rho_mat * Ev;
    
    if (dump_id == total_num_dumps)
        final_rho_and = current_rho_mat_and;
    end
    
    if dump_id > 1
        local_entropies(dump_id) = -trace(current_rho_mat * logm(current_rho_mat));
        local_entropies_and(dump_id) = -trace(current_rho_mat_and * logm(current_rho_mat_and));
    end
    
    n_part = zeros(1,Nc);
    n_part_and = zeros(1,Nc);
    for k = 1:Ns
        n_part = n_part + real(current_rho_mat(k,k)) * ( dec2bin(idtox(k),Nc) == '1' );
        n_part_and = n_part_and + real(current_rho_mat_and(k,k)) * ( dec2bin(idtox(k),Nc) == '1' );
    end
    
    local_imbalances(dump_id) = sum(n_part(1:2:Nc-1) - n_part(2:2:Nc)) / sum(n_part);
    local_imbalances_and(dump_id) = sum(n_part_and(1:2:Nc-1) - n_part_and(2:2:Nc)) / sum(n_part_and);
    
    stationary_ipr = 0.0;
    stationary_ipr_and = 0.0;
    for s_id = 1:Ns
        stationary_ipr = stationary_ipr + abs(current_rho_mat(s_id,s_id))^2;
        stationary_ipr_and = stationary_ipr_and + abs(current_rho_mat_and(s_id,s_id))^2;
    end
    
    local_iprs(dump_id) = stationary_ipr;
    local_iprs_and(dump_id) = stationary_ipr_and;
end

total_num_dumps = size(local_imbalances, 1);

file_name = sprintf('%sentropy_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', local_entropies(dump_id));
end
fclose(file_id);

file_name = sprintf('%sentropy_and_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', local_entropies_and(dump_id));
end
fclose(file_id);

file_name = sprintf('%simbalance_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', local_imbalances(dump_id));
end
fclose(file_id);

file_name = sprintf('%simbalance_and_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', local_imbalances_and(dump_id));
end
fclose(file_id);

file_name = sprintf('%sipr_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', real(local_iprs(dump_id)));
end
fclose(file_id);

file_name = sprintf('%sipr_and_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:total_num_dumps
    fprintf(file_id, '%0.18e\n', real(local_iprs_and(dump_id)));
end
fclose(file_id);

if(save_type >= 1)
    file_name = sprintf('%sfinal_rho_diag_in_direct_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for s_id = 1:Ns
        fprintf(file_id, '%0.18e\n', final_rho(s_id, s_id));
    end
    fclose(file_id);
end

if(save_type >= 2)
	file_name = sprintf('%sfinal_rho_in_direct_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
	for s_id_1 = 1:Ns
		for s_id_2 = 1:Ns
			fprintf(file_id, '%0.18e %0.18e\n', real(final_rho(s_id_1, s_id_2)), imag(final_rho(s_id_1, s_id_2)) );
		end  
	end
	fclose(file_id);
end

if(save_type >= 1)
    file_name = sprintf('%sfinal_rho_diag_in_stationary_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for s_id = 1:Ns
        fprintf(file_id, '%0.18e\n', final_rho_and(s_id, s_id));
    end
    fclose(file_id);
end

if(save_type >= 2)
	file_name = sprintf('%sfinal_rho_in_stationary_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
	for s_id_1 = 1:Ns
		for s_id_2 = 1:Ns
			fprintf(file_id, '%0.18e %0.18e\n', real(final_rho_and(s_id_1, s_id_2)), imag(final_rho_and(s_id_1, s_id_2)) );
		end  
	end
	fclose(file_id);
end

[zero_evec,zero_eval] = eigs(lndbldn_sprs, 1, 'sm');
stationary_rho_array = zero_evec;

st_rho = zeros(Ns, Ns);
for s_id=1:Ns
    st_rho(:,s_id) = stationary_rho_array(1+(s_id-1)*(Ns) : s_id*(Ns));
end
st_rho = st_rho/trace(st_rho);

st_rho_and = Ev'*st_rho*Ev;

diag_diff = zeros(Ns, 1);
diag_diff_rel = zeros(Ns, 1);
diag_st_rho = diag(st_rho);
diag_final_rho = diag(final_rho);
for s_id = 1:Ns
    diag_diff(s_id) = abs(abs(diag_st_rho(s_id)) - abs(diag_final_rho(s_id)));
    diag_diff_rel(s_id) = diag_diff(s_id) / abs(diag_st_rho(s_id));
end

stationary_entropy = -trace(st_rho * logm(st_rho)); % entropy is direct basis
stationary_entropy_and = -trace(st_rho_and * logm(st_rho_and)); % entropy in stationary basis

n_part = zeros(1,Nc);
n_part_and = zeros(1,Nc);
for s_id=1:Ns
    n_part = n_part + real(st_rho(s_id,s_id))*(dec2bin(idtox(s_id),Nc)=='1');
    n_part_and = n_part_and + real(st_rho_and(s_id,s_id))*(dec2bin(idtox(s_id),Nc)=='1');
end

stationary_imbalance = sum(n_part(1:2:Nc-1)-n_part(2:2:Nc)) / sum(n_part);
stationary_imbalance_and = sum(n_part_and(1:2:Nc-1)-n_part_and(2:2:Nc)) / sum(n_part_and);

stationary_ipr = 0.0;
stationary_ipr_and = 0.0;
for s_id = 1:Ns
    stationary_ipr = stationary_ipr + abs(st_rho(s_id,s_id))^2;
    stationary_ipr_and = stationary_ipr_and + abs(st_rho_and(s_id,s_id))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Checking Discrepancy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

discrepancy = -sqrt(-1) * (H*st_rho - st_rho*H);

if(diss_type == 1)
    for diss_id = 1:Nc-1
        diss = zeros(Ns);
        for s_id_1 = 1:Ns
            diss(s_id_1, s_id_1) = bitget(idtox(s_id_1), diss_id) - bitget(idtox(s_id_1), diss_id+1);
            for s_id_2 = 1:Ns
                if(is_adjacent(idtox(s_id_1), idtox(s_id_2)))
                    hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(s_id_1), idtox(s_id_2)),Nc) == '1' );
                    if(hopping_ids(2) == diss_id)
                        if(bitget(idtox(s_id_1), diss_id))
                            diss(s_id_2,s_id_1) = exp(sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = -exp(-sqrt(-1)*diss_phase);
                        else
                            diss(s_id_2,s_id_1) = -exp(-sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = exp(sqrt(-1)*diss_phase);
                        end
                    end
                end
            end
        end
        
        discrepancy = discrepancy + g * 0.5 * ...
            (2.0 * diss * st_rho * diss' - ...
            st_rho * diss' * diss - ...
            diss' * diss * st_rho);
    end
    
    if (periodic_bc == 1)
        
        diss = zeros(Ns);
        for s_id_1 = 1:Ns
            diss(s_id_1, s_id_1) = bitget(idtox(s_id_1), Nc) - bitget(idtox(s_id_1), 1);
            for s_id_2 = 1:Ns
                if(is_adjacent(idtox(s_id_1), idtox(s_id_2)))
                    hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(s_id_1), idtox(s_id_2)),Nc) == '1' );
                    if(hopping_ids(2) == Nc)
                        if(bitget(idtox(s_id_1), Nc))
                            diss(s_id_2,s_id_1) = exp(sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = -exp(-sqrt(-1)*diss_phase);
                        else
                            diss(s_id_2,s_id_1) = -exp(-sqrt(-1)*diss_phase);
                            diss(s_id_1,s_id_2) = exp(sqrt(-1)*diss_phase);
                        end
                    end
                end
            end
        end
        
        discrepancy = discrepancy + g * 0.5 * ...
            (2.0 * diss * st_rho * diss' - ...
            st_rho * diss' * diss - ...
            diss' * diss * st_rho);
        
    end
end

if(diss_type == 0)
    for diss_id = 1:Nc
        diss = zeros(Ns);
        for s_id_1 = 1:Ns
            diss(s_id_1, s_id_1) = bitget(idtox(s_id_1), diss_id);
        end
        
        discrepancy = discrepancy + g * 0.5 * ...
            (2.0 * diss * st_rho * diss' - ...
            st_rho * diss' * diss - ...
            diss' * diss * st_rho);
    end
end

error = max(max(abs(discrepancy)));

file_name = sprintf('%sstationary_info_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
fprintf(file_id, '%0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e\n', error, abs(zero_eval), max(diag_diff), max(diag_diff_rel), stationary_entropy, stationary_entropy_and, stationary_imbalance, stationary_imbalance_and, stationary_ipr, stationary_ipr_and);
fclose(file_id);

file_name = sprintf('%sstationary_rho_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
st_rho_abs_diag = abs(diag(st_rho));
for s_id = 1:Ns
    fprintf(file_id, '%0.18e\n', st_rho_abs_diag(s_id));
end
fclose(file_id);

