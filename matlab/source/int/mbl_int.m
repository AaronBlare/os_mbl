clear all;

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
seed_start          = input_data(10);
seed_num            = input_data(11);
is_int              = input_data(12);
int_ist             = input_data(13);
int_isi             = input_data(14);
int_dt              = input_data(15);
int_db              = input_data(16);
int_de              = input_data(17);
int_dn              = input_data(18);
is_zev              = input_data(19);
is_zev_check        = input_data(20);
is_eg_stat          = input_data(21);
is_save_vec         = input_data(22);
is_save_mtx         = input_data(23);
fs_type             = input_data(24);

if fs_type == 1
    data_path = '../../../data/matlab/';
elseif fs_type == 0
    data_path = '';
else
    error('Error: wrong dump_type');
end

Np = Nc/2;
Ns = nchoosek(Nc, Np);
num = precalc_states (Nc, Np, periodic_bc);

for seed = seed_start : seed_start + (seed_num - 1)
    
    fn_suffix = sprintf('Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disorder
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Hd = zeros(Ns);
    rng(seed)
    E = 2 * rand(Nc,1) - 1;
    
    if (energy_type == 1)
        E = E - sum(E)/Nc;
    end
    
    for s_id = 1:Ns
        Hd(s_id,s_id) = ( dec2bin(idtox(s_id), Nc) == '1' ) * E;
    end
    
    file_name = sprintf('%srandom_energies_%s', data_path, fn_suffix);
    file_id = fopen(file_name, 'w');
    for c_id = 1:Nc
        fprintf(file_id, '%0.18e\n', E(c_id));
    end
    fclose(file_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interaction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hopping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Hh = zeros(Ns);
    for s_id_1 = 1:Ns
        for s_id_2 = 1:Ns
            Hh(s_id_1,s_id_2) = is_adjacent(idtox(s_id_1), idtox(s_id_2));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Total Hamiltonian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H= -J*Hh + U*Hi + 2.0*W*Hd;
    
    if(is_save_mtx == 1)
        file_name = sprintf('%shamiltonian_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for s_id_1 = 1:Ns
            for s_id_2 = 1:Ns
                fprintf(file_id, '%0.18e\n', H(s_id_1, s_id_2));
            end
        end
        fclose(file_id);
    end
    
    [Ev,Eg] = eig(H);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Saving
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(is_save_vec == 1)
        file_name = sprintf('%seg_hamiltonian_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for s_id = 1:Ns
            fprintf(file_id, '%0.18e\n', Eg(s_id, s_id));
        end
        fclose(file_id);
    end
    
    if(is_save_mtx == 1)
        file_name = sprintf('%sev_hamiltonian_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for s_id_1 = 1:Ns
            for s_id_2 = 1:Ns
                fprintf(file_id, '%0.18e\n', Ev(s_id_1, s_id_2));
            end
        end
        fclose(file_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creating matrix of right part lndbldn*Rho=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lndbldn = zeros((Ns)^2);
    lndbldn = -sqrt(-1) * (kron(eye(Ns),H) - kron(transpose(H),eye(Ns)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dissipators
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
          
    elseif (diss_type == 2)
        
        diss_id = Nc/2;
        
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sparse Lindbladian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lndbldn_sprs = sparse(lndbldn);
    lndbldn = 0;
    
    if (is_int == 1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Init Integration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fn_suffix = sprintf('int_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d)_iist(%d)_iisi(%d)_idt(%d).txt', ...
            Nc, ...
            diss_type, ...
            diss_phase, ...
            energy_type, ...
            periodic_bc, ...
            W, ...
            U, ...
            J, ...
            g, ...
            seed, ...
            int_ist, ...
            int_isi, ...
            int_dt);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Integration times
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dump_times = zeros(int_dn + 1, 1);
        if (int_dt == 0)
            dump_shift = (int_de - int_db) / int_dn;
            for dump_id = 1:int_dn+1
                dump_times(dump_id) = int_db + (dump_id - 1) * dump_shift;
            end
        elseif (int_dt == 1)
            begin_decade = log10(int_db);
            end_decade = log10(int_de);
            num_decades = end_decade - begin_decade;
            num_decade_dumps = int_dn / num_decades;
            for dump_id = 1:int_dn+1
                dump_times(dump_id) = power(10, begin_decade) * power(10, (1.0/num_decade_dumps) * (dump_id - 1));
            end
        else
            error('Error: wrong dump_type');
        end
        
        dump_times = vertcat(0, dump_times);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Integration Init Conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        start_rho = zeros(Ns * Ns, 1);
        if int_ist == 0
            start_rho((int_isi - 1) * (Ns) + int_isi) = 1.0;
        elseif int_ist == 1
            for s_id = 1:Ns
                start_rho((s_id - 1) * (Ns) + s_id) = 1.0 / Ns;
            end
        else
            error('Error: wrong init_state_type');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Integration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tic
        [times, rho_dumps] = ode45(@(times, rho_dumps) right_part(times, rho_dumps, lndbldn_sprs), dump_times, start_rho);
        toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Integration Calc Characteristics Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        total_num_dumps = size(times, 1);
        
        local_entropies         = zeros(total_num_dumps, 1);
        local_entropies_st      = zeros(total_num_dumps, 1);
        local_imbalances        = zeros(total_num_dumps, 1);
        local_imbalances_and    = zeros(total_num_dumps, 1);
        local_iprs              = zeros(total_num_dumps, 1);
        local_iprs_st           = zeros(total_num_dumps, 1);
        
        int_rho = zeros(Ns, Ns);
        int_rho_st = zeros(Ns, Ns);
        
        for dump_id = 1:total_num_dumps
            
            current_rho_vec = rho_dumps(dump_id, :);
            
            current_rho_mat = zeros(Ns,Ns);
            for s_id_1 = 1:Ns
                for s_id_2 = 1:Ns
                    current_rho_mat(s_id_1, s_id_2) = current_rho_vec((s_id_1-1)*(Ns) + s_id_2);
                end
            end
            
            if (dump_id == total_num_dumps)
                int_rho = current_rho_mat;
            end
            
            current_rho_mat_and = Ev' * current_rho_mat * Ev;
            
            if (dump_id == total_num_dumps)
                int_rho_st = current_rho_mat_and;
            end
            
            if dump_id > 1
                local_entropies(dump_id) = -trace(current_rho_mat * logm(current_rho_mat));
                local_entropies_st(dump_id) = -trace(current_rho_mat_and * logm(current_rho_mat_and));
            end
            
            n_part = zeros(1,Nc);
            n_part_st = zeros(1,Nc);
            for k = 1:Ns
                n_part = n_part + real(current_rho_mat(k,k)) * ( dec2bin(idtox(k),Nc) == '1' );
                n_part_st = n_part_st + real(current_rho_mat_and(k,k)) * ( dec2bin(idtox(k),Nc) == '1' );
            end
            
            local_imbalances(dump_id) = sum(n_part(1:2:Nc-1) - n_part(2:2:Nc)) / sum(n_part);
            local_imbalances_and(dump_id) = sum(n_part_st(1:2:Nc-1) - n_part_st(2:2:Nc)) / sum(n_part_st);
            
            zev_ipr = 0.0;
            zev_ipr_st = 0.0;
            for s_id = 1:Ns
                zev_ipr = zev_ipr + abs(current_rho_mat(s_id,s_id))^2;
                zev_ipr_st = zev_ipr_st + abs(current_rho_mat_and(s_id,s_id))^2;
            end
            
            local_iprs(dump_id) = zev_ipr;
            local_iprs_st(dump_id) = zev_ipr_st;
        end
        
        if (is_eg_stat == 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eigen Values Statistics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            eval_rho = eig(int_rho);
            
            file_name = sprintf('%srho_evals_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id_1 = 1:Ns
                fprintf(file_id, '%0.18e\n', eval_rho(s_id_1));
            end
            fclose(file_id);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Integration Dump Characteristics Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        file_name = sprintf('%stimes_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', times(dump_id));
        end
        fclose(file_id);
        
        total_num_dumps = size(local_imbalances, 1);
        
        file_name = sprintf('%sentropy_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', local_entropies(dump_id));
        end
        fclose(file_id);
        
        file_name = sprintf('%sentropy_st_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', local_entropies_st(dump_id));
        end
        fclose(file_id);
        
        file_name = sprintf('%simbalance_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', local_imbalances(dump_id));
        end
        fclose(file_id);
        
        file_name = sprintf('%simbalance_st_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', local_imbalances_and(dump_id));
        end
        fclose(file_id);
        
        file_name = sprintf('%sipr_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', real(local_iprs(dump_id)));
        end
        fclose(file_id);
        
        file_name = sprintf('%sipr_st_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        for dump_id = 1:total_num_dumps
            fprintf(file_id, '%0.18e\n', real(local_iprs_st(dump_id)));
        end
        fclose(file_id);
        
        if(is_save_vec == 1)
            file_name = sprintf('%srho_diag_d_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id = 1:Ns
                fprintf(file_id, '%0.18e\n', int_rho(s_id, s_id));
            end
            fclose(file_id);
        end
        
        if(is_save_mtx >= 1)
            file_name = sprintf('%srho_d_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id_1 = 1:Ns
                for s_id_2 = 1:Ns
                    fprintf(file_id, '%0.18e %0.18e\n', real(int_rho(s_id_1, s_id_2)), imag(int_rho(s_id_1, s_id_2)) );
                end
            end
            fclose(file_id);
        end
        
        if(is_save_vec == 1)
            file_name = sprintf('%srho_diag_st_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id = 1:Ns
                fprintf(file_id, '%0.18e\n', int_rho_st(s_id, s_id));
            end
            fclose(file_id);
        end
        
        if(is_save_mtx >= 1)
            file_name = sprintf('%srho_st_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id_1 = 1:Ns
                for s_id_2 = 1:Ns
                    fprintf(file_id, '%0.18e %0.18e\n', real(int_rho_st(s_id_1, s_id_2)), imag(int_rho_st(s_id_1, s_id_2)) );
                end
            end
            fclose(file_id);
        end
        
    end
    
    if (is_zev == 1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Zero Eigen Vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fn_suffix = sprintf('zev_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
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
        
        tic
        
        [zero_evec,zero_eval] = eigs(lndbldn_sprs, 1, 'sm');
        stationary_rho_array = zero_evec;
        
        toc
        
        zev_rho = zeros(Ns, Ns);
        for s_id=1:Ns
            zev_rho(:,s_id) = stationary_rho_array(1+(s_id-1)*(Ns) : s_id*(Ns));
        end
        zev_rho = zev_rho/trace(zev_rho);
        
        zev_rho_st = Ev'*zev_rho*Ev;
        
        if (is_eg_stat == 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eigen Values Statistics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            eval_rho = eig(zev_rho);
            
            file_name = sprintf('%srho_evals_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id_1 = 1:Ns
                fprintf(file_id, '%0.18e\n', eval_rho(s_id_1));
            end
            fclose(file_id);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Zero Eigen Vector Calc Characteristics Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        zev_entropy = -trace(zev_rho * logm(zev_rho));
        zev_entropy_st = -trace(zev_rho_st * logm(zev_rho_st));
        
        n_part = zeros(1,Nc);
        n_part_st = zeros(1,Nc);
        for s_id=1:Ns
            n_part = n_part + real(zev_rho(s_id,s_id))*(dec2bin(idtox(s_id),Nc)=='1');
            n_part_st = n_part_st + real(zev_rho_st(s_id,s_id))*(dec2bin(idtox(s_id),Nc)=='1');
        end
        
        zev_imbalance = sum(n_part(1:2:Nc-1)-n_part(2:2:Nc)) / sum(n_part);
        zev_imbalance_st = sum(n_part_st(1:2:Nc-1)-n_part_st(2:2:Nc)) / sum(n_part_st);
        
        zev_ipr = 0.0;
        zev_ipr_st = 0.0;
        for s_id = 1:Ns
            zev_ipr = zev_ipr + abs(zev_rho(s_id,s_id))^2;
            zev_ipr_st = zev_ipr_st + abs(zev_rho_st(s_id,s_id))^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Zero Eigen Vector Dump Characteristics Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        file_name = sprintf('%sentropy_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', zev_entropy);
        fclose(file_id);
        
        file_name = sprintf('%sentropy_st_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', zev_entropy_st);
        fclose(file_id);
        
        file_name = sprintf('%simbalance_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', zev_imbalance);
        fclose(file_id);
        
        file_name = sprintf('%simbalance_st_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', zev_imbalance_st);
        fclose(file_id);
        
        file_name = sprintf('%sipr_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', zev_ipr);
        fclose(file_id);
        
        file_name = sprintf('%sipr_st_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', zev_ipr_st);
        fclose(file_id);
        
        if(is_save_vec == 1)
            file_name = sprintf('%srho_diag_d_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id = 1:Ns
                fprintf(file_id, '%0.18e\n', zev_rho(s_id, s_id));
            end
            fclose(file_id);
        end
        
        if(is_save_mtx >= 1)
            file_name = sprintf('%srho_d_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id_1 = 1:Ns
                for s_id_2 = 1:Ns
                    fprintf(file_id, '%0.18e %0.18e\n', real(zev_rho(s_id_1, s_id_2)), imag(zev_rho(s_id_1, s_id_2)) );
                end
            end
            fclose(file_id);
        end
        
        if(is_save_vec == 1)
            file_name = sprintf('%srho_diag_st_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id = 1:Ns
                fprintf(file_id, '%0.18e\n', zev_rho_st(s_id, s_id));
            end
            fclose(file_id);
        end
        
        if(is_save_mtx >= 1)
            file_name = sprintf('%srho_st_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            for s_id_1 = 1:Ns
                for s_id_2 = 1:Ns
                    fprintf(file_id, '%0.18e %0.18e\n', real(zev_rho_st(s_id_1, s_id_2)), imag(zev_rho_st(s_id_1, s_id_2)) );
                end
            end
            fclose(file_id);
        end
        
        if (is_zev_check == 1)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Checking Discrepancy
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            discrepancy = -sqrt(-1) * (H*zev_rho - zev_rho*H);
            
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
                        (2.0 * diss * zev_rho * diss' - ...
                        zev_rho * diss' * diss - ...
                        diss' * diss * zev_rho);
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
                        (2.0 * diss * zev_rho * diss' - ...
                        zev_rho * diss' * diss - ...
                        diss' * diss * zev_rho);
                    
                end
            end
            
            if(diss_type == 2)
                diss_id = Nc/2;
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
                    (2.0 * diss * zev_rho * diss' - ...
                    zev_rho * diss' * diss - ...
                    diss' * diss * zev_rho);
            end
            
            if(diss_type == 0)
                for diss_id = 1:Nc
                    diss = zeros(Ns);
                    for s_id_1 = 1:Ns
                        diss(s_id_1, s_id_1) = bitget(idtox(s_id_1), diss_id);
                    end
                    
                    discrepancy = discrepancy + g * 0.5 * ...
                        (2.0 * diss * zev_rho * diss' - ...
                        zev_rho * diss' * diss - ...
                        diss' * diss * zev_rho);
                end
            end
            
            error = max(max(abs(discrepancy)));
            
            file_name = sprintf('%serror_%s', data_path, fn_suffix);
            file_id = fopen(file_name, 'w');
            fprintf(file_id, '%0.18e\n', error);
            fprintf(file_id, '%0.18e\n', abs(zero_eval));
            fclose(file_id);
        end
    end
    
    if (is_int && is_zev)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Intergation vs Zero Eigen Vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fn_suffix = sprintf('int_zev_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt', ...
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
        
        diag_diff_abs = zeros(Ns, 1);
        diag_diff_rel = zeros(Ns, 1);
        diag_int_rho = diag(int_rho);
        diag_zev_rho = diag(zev_rho);
        for s_id = 1:Ns
            diag_diff_abs(s_id) = abs(abs(diag_zev_rho(s_id)) - abs(diag_int_rho(s_id)));
            diag_diff_rel(s_id) = diag_diff_abs(s_id) / abs(diag_zev_rho(s_id));
        end
        
        
        file_name = sprintf('%sdiff_%s', data_path, fn_suffix);
        file_id = fopen(file_name, 'w');
        fprintf(file_id, '%0.18e\n', max(diag_diff_abs));
        fprintf(file_id, '%0.18e\n', max(diag_diff_rel));
        fclose(file_id);
    end
    
end
