clear all;

tic

filename = 'config.txt';
input_data = importdata(filename);

Nc                  = input_data(1);
dissipator_type     = input_data(2);
alpha               = input_data(3);
energy_type         = input_data(4);
W                   = input_data(5);
U                   = input_data(6);
J                   = input_data(7);
g                   = input_data(8);
seed                = input_data(9);
init_state_type     = input_data(10);
init_state_id       = input_data(11);
dump_type           = input_data(12);
begin_dump          = input_data(13);
end_dump            = input_data(14);
num_dumps           = input_data(15);
save_type           = input_data(16);
dump_fs           = input_data(17);

if dump_fs == 0
    data_path = '../../../data/int/matlab/';
elseif dump_fs == 1
    data_path = '';
else
   error('Error: wrong dump_type');
end

Np = Nc/2;
Ns = nchoosek(Nc, Np);
num = precalc_states (Nc, Np);

file_name_suffix = sprintf('Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ist(%d)_iss(%d)_dump_type(%d)_seed(%d).txt', ...
    Nc, ...
    dissipator_type, ...
    alpha, ...
    energy_type, ...
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

for k = 1:Ns
    Hd(k,k) = ( dec2bin(idtox(k), Nc) == '1' ) * E;
end

file_name = sprintf('%srandom_energies_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
for dump_id = 1:Nc
    fprintf(file_id, '%0.18e\n', E(dump_id));
end
fclose(file_id);

%%%%%%%%%%%%%%%%%%%%%
% Interaction
%%%%%%%%%%%%%%%%%%%%%
Hi = zeros(Ns);
for k = 1:Ns
    Hi(k,k) = sum(dec2bin(bitand(idtox(k), bitshift(idtox(k), 1))) == '1');
end

%%%%%%%%%%%%%%%%%%%%%
% Hopping
%%%%%%%%%%%%%%%%%%%%%
Hh = zeros(Ns);
for k = 1:Ns
    for kk = 1:Ns
        Hh(k,kk) = is_adjacent(idtox(k), idtox(kk));
    end
end

%%%%%%%%%%%%
% Total
%%%%%%%%%%%%
H= -J*Hh + U*Hi + 2.0*W*Hd;

if(save_type >= 1)
    file_name = sprintf('%shamiltonian_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for state_id_1 = 1:Ns
        for state_id_2 = 1:Ns
            fprintf(file_id, '%0.18e\n', H(state_id_1, state_id_2));
        end
    end
    fclose(file_id);
end

[Ev,Eg] = eig(H); % Anderson modes

if(save_type >= 1)
    file_name = sprintf('%seg_hamiltonian_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for state_id = 1:Ns
        fprintf(file_id, '%0.18e\n', Eg(state_id, state_id));
    end
    fclose(file_id);
end

if(save_type >= 2)
    file_name = sprintf('%sev_hamiltonian_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for state_id_1 = 1:Ns
        for state_id_2 = 1:Ns
            fprintf(file_id, '%0.18e\n', Ev(state_id_1, state_id_2));
        end
    end
    fclose(file_id);
end

%%%%%%%%%%%%%%%%%%%
% Creating matrix of right part P*Rho=0
%%%%%%%%%%%%%%%%%%%
P = zeros((Ns)^2);
P=-sqrt(-1)*(kron(eye(Ns),H)-kron(transpose(H),eye(Ns)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (dissipator_type == 1)
    
    for k=1:Nc-1
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k)-bitget(idtox(kk),k+1);
            for kkk=1:Ns
                if(is_adjacent(idtox(kk), idtox(kkk)))
                    hop=1+Nc-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),Nc)=='1');
                    if(hop(2)==k)
                        if((bitget(idtox(kk),k)))
                            A(kkk,kk) = 1.0 * exp(sqrt(-1)*alpha);
                            A(kk,kkk) = -1.0 * exp(-sqrt(-1)*alpha);
                        else
                            A(kkk,kk) = -1.0 * exp(-sqrt(-1)*alpha);
                            A(kk,kkk) = 1.0 * exp(sqrt(-1)*alpha);
                        end
                    end
                end
            end
        end
        
        P=P+...
            g/2*(2*kron(eye(Ns),A)*kron(transpose(A'),eye(Ns))-...
            kron(transpose(A'*A),eye(Ns))-kron(eye(Ns),A'*A));
        
    end
elseif(dissipator_type == 0)
    
    for k=1:Nc
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k);
        end
        
        P=P+...
            g/2*(2*kron(eye(Ns),A)*kron(transpose(A'),eye(Ns))-...
            kron(transpose(A'*A),eye(Ns))-kron(eye(Ns),A'*A));
        
    end
else
    error('Error: wrong dissipator_type');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

Ps = sparse(P);
P=0;

start_rho = zeros(Ns * Ns, 1);
if init_state_type == 0
    start_rho((init_state_id - 1) * (Ns) + init_state_id) = 1.0;
elseif init_state_type == 1
    for i = 1:Ns
        start_rho((i - 1) * (Ns) + i) = 1.0 / Ns;
    end
else
    error('Error: wrong init_state_type');
end

tic
[times, rho_dumps] = ode45(@(times, rho_dumps) right_part(times, rho_dumps, Ps), dump_times, start_rho);
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
    for i = 1:Ns
        for j = 1:Ns
            current_rho_mat(i, j) = current_rho_vec((i-1)*(Ns) + j);
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
    for i = 1:Ns
        stationary_ipr = stationary_ipr + abs(current_rho_mat(i,i))^2;
        stationary_ipr_and = stationary_ipr_and + abs(current_rho_mat_and(i,i))^2;
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
    for state_id = 1:Ns
        fprintf(file_id, '%0.18e\n', final_rho(state_id, state_id));
    end
    fclose(file_id);
end

if(save_type >= 2)
	file_name = sprintf('%sfinal_rho_in_direct_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
	for state_id_1 = 1:Ns
		for state_id_2 = 1:Ns
			fprintf(file_id, '%0.18e %0.18e\n', real(final_rho(state_id_1, state_id_2)), imag(final_rho(state_id_1, state_id_2)) );
		end  
	end
	fclose(file_id);
end

if(save_type >= 1)
    file_name = sprintf('%sfinal_rho_diag_in_stationary_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for state_id = 1:Ns
        fprintf(file_id, '%0.18e\n', final_rho_and(state_id, state_id));
    end
    fclose(file_id);
end

if(save_type >= 2)
	file_name = sprintf('%sfinal_rho_in_stationary_basis_%s', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
	for state_id_1 = 1:Ns
		for state_id_2 = 1:Ns
			fprintf(file_id, '%0.18e %0.18e\n', real(final_rho_and(state_id_1, state_id_2)), imag(final_rho_and(state_id_1, state_id_2)) );
		end  
	end
	fclose(file_id);
end

[zero_evec,zero_eval] = eigs(Ps, 1, 'sm');
stationary_rho_array = zero_evec;

st_rho = zeros(Ns, Ns);
for i=1:Ns
    st_rho(:,i) = stationary_rho_array(1+(i-1)*(Ns) : i*(Ns));
end
st_rho = st_rho/trace(st_rho);

st_rho_and = Ev'*st_rho*Ev;

diag_diff = zeros(Ns, 1);
diag_diff_rel = zeros(Ns, 1);
diag_st_rho = diag(st_rho);
diag_final_rho = diag(final_rho);
for i = 1:Ns
    diag_diff(i) = abs(abs(diag_st_rho(i)) - abs(diag_final_rho(i)));
    diag_diff_rel(i) = diag_diff(i) / abs(diag_st_rho(i));
end

stationary_entropy = -trace(st_rho * logm(st_rho)); % entropy is direct basis
stationary_entropy_and = -trace(st_rho_and * logm(st_rho_and)); % entropy in stationary basis

n_part = zeros(1,Nc);
n_part_and = zeros(1,Nc);
for k=1:Ns
    n_part = n_part + real(st_rho(k,k))*(dec2bin(idtox(k),Nc)=='1');
    n_part_and = n_part_and + real(st_rho_and(k,k))*(dec2bin(idtox(k),Nc)=='1');
end

stationary_imbalance = sum(n_part(1:2:Nc-1)-n_part(2:2:Nc)) / sum(n_part);
stationary_imbalance_and = sum(n_part_and(1:2:Nc-1)-n_part_and(2:2:Nc)) / sum(n_part_and);

stationary_ipr = 0.0;
stationary_ipr_and = 0.0;
for i = 1:Ns
    stationary_ipr = stationary_ipr + abs(st_rho(i,i))^2;
    stationary_ipr_and = stationary_ipr_and + abs(st_rho_and(i,i))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=-sqrt(-1)*(H*st_rho-st_rho*H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Dissipators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(dissipator_type == 1)
    for k = 1:Nc-1
        A = zeros(Ns);
        for kk = 1:Ns
            A(kk,kk) = bitget(idtox(kk),k) - bitget(idtox(kk),k+1);
            for kkk = 1:Ns
                if(is_adjacent(idtox(kk), idtox(kkk)))
                    hop=1+Nc-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),Nc)=='1');
                    if(hop(2)==k)
                        
                        if((~bitget(idtox(kk),k))&&(bitget(idtox(kkk),k+1)))
                            A(kkk,kk)=1;
                            A(kk,kkk)=-1;
                        else
                            A(kkk,kk)=-1;
                            A(kk,kkk)=1;
                        end
                    end
                end
            end
        end
        
        dy=dy+g/2*(2*A*st_rho*A'-st_rho*A'*A-A'*A*st_rho);
    end
elseif(dissipator_type == 0)
    
    for k = 1:Nc
        A = zeros(Ns);
        for kk = 1:Ns
            A(kk,kk) = bitget(idtox(kk),k);
        end
        
        dy=dy+g/2*(2*A*st_rho*A' - st_rho*A'*A - A'*A*st_rho);
    end
else
    error('Error: wrong dissipator_type');
end

err = max(max(abs(dy)));

file_name = sprintf('%sstationary_info_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
fprintf(file_id, '%0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e %0.18e\n', err, abs(zero_eval), max(diag_diff), max(diag_diff_rel), stationary_entropy, stationary_entropy_and, stationary_imbalance, stationary_imbalance_and, stationary_ipr, stationary_ipr_and);
fclose(file_id);

file_name = sprintf('%sstationary_rho_%s', data_path, file_name_suffix);
file_id = fopen(file_name, 'w');
st_rho_abs_diag = abs(diag(st_rho));
for i = 1:Ns
    fprintf(file_id, '%0.18e\n', st_rho_abs_diag(i));
end
fclose(file_id);

