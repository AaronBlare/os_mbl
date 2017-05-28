clear all;

filename = 'config.txt';
input_data = importdata(filename);

Nc                  = input_data(1);
dissipator_type     = input_data(2);
alpha               = input_data(3);		
energy_type         = input_data(4);
border_conditions   = input_data(5);
W                   = input_data(6);
U                   = input_data(7);
J                   = input_data(8);
g                   = input_data(9);
seed_start          = input_data(10);
seed_num            = input_data(11);
save_type           = input_data(12);
save_super_operator = input_data(13);
dump_type           = input_data(14);

if dump_type == 1
    data_path = '../../../data/zev/matlab/';
elseif dump_type == 0
    data_path = '';
else
   error('Error: wrong dump_type');
end

Np = Nc/2;
Ns = nchoosek(Nc,Np);
num_states = precalc_states(Nc, Np, border_conditions);

for seed = seed_start : seed_start + (seed_num - 1)

    tic
    
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
    seed_start, ...
    seed_num, ...
    seed);
    
    H  = zeros(Ns); % Halimtonian
    Hd = zeros(Ns); % Disorder
    Hi = zeros(Ns); % Interaction
    Hh = zeros(Ns); % Hopping
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             Disorder
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rng(seed)
    E = 2 * rand(Nc,1) - 1;
    if energy_type == 1
        E = E - sum(E) / Nc;
    end
    
    for state_id_1 = 1:Ns
        Hd(state_id_1, state_id_1) = ( dec2bin(idtox(state_id_1), Nc) == '1' ) * E;
    end
    
    file_name = sprintf('%srandom_energies_%s.txt', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    for dump_id = 1:Nc
        fprintf(file_id, '%0.18e\n', E(dump_id));
    end
    fclose(file_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            Interaction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for state_id_1 = 1:Ns
        shifted = bitshift(idtox(state_id_1), 1);
        
        if (border_conditions == 1)
            if (idtox(state_id_1) >= 2^(Nc-1))
                shifted = shifted + 1;
            end
        end
        
        Hi(state_id_1, state_id_1) = sum( dec2bin(bitand(idtox(state_id_1), shifted)) == '1' );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                              Hopping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for state_id_1 = 1:Ns
        for state_id_2 = 1:Ns
            Hh(state_id_1, state_id_2) = is_adjacent(idtox(state_id_1), idtox(state_id_2));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            Hamiltonian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    H = -J*Hh + U*Hi + 2.0*W*Hd;
    
    if(save_type >= 2)
        file_name = sprintf('%shamiltonian_%s.txt',data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e\n', H(state_id_1, state_id_2));
            end
        end
        fclose(file_id);
    end
    
    Hd = 0;
    Hi = 0;
    Hh = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     Hamiltonian eigen problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [evecs_H, evals_H] = eig(H); % Anderson modes
    evals_H = diag(evals_H);
    
    if(save_type >= 1)
        file_name = sprintf('%shamiltonian_evals_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            fprintf(file_id, '%0.18e\n', evals_H(state_id_1));
        end
        fclose(file_id);
    end
    
    if(save_type >= 2)
        file_name = sprintf('%shamiltonian_evecs_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e\n', evecs_H(state_id_1, state_id_2));
            end
        end
        fclose(file_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                Creating supermatrix of right part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    super_rp_matrix = -sqrt(-1) * ( kron(eye(Ns), H) - kron(transpose(H), eye(Ns)) );
    
    if(save_super_operator >= 0)
        file_name = sprintf('%slindbladian_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1 : Ns * Ns
            for state_id_2 = 1 : Ns * Ns
                fprintf(file_id, '%0.18e %0.18e\n', real(super_rp_matrix(state_id_1, state_id_2)), imag(super_rp_matrix(state_id_1, state_id_2)));
            end
        end
        fclose(file_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             Dissipator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(dissipator_type == 1)
        for dissipator_id = 1:Nc-1
            dissipator = zeros(Ns);
            for state_id_1 = 1:Ns
                dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), dissipator_id) - bitget(idtox(state_id_1), dissipator_id+1);
                for state_id_2 = 1:Ns
                    if(is_adjacent(idtox(state_id_1), idtox(state_id_2)))
                        hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(state_id_1), idtox(state_id_2)),Nc) == '1' );
                        if(hopping_ids(2) == dissipator_id)
                            if(bitget(idtox(state_id_1), dissipator_id))
                                dissipator(state_id_2,state_id_1) = exp(sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = -exp(-sqrt(-1)*alpha);
                            else
                                dissipator(state_id_2,state_id_1) = -exp(-sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = exp(sqrt(-1)*alpha);
                            end
                        end
                    end
                end
            end
            
            if(save_super_operator >= 0)
                file_name = sprintf('%sdissipator_%d_%s.txt', data_path, dissipator_id, file_name_suffix);
                file_id = fopen(file_name, 'w');
                for state_id_1 = 1 : Ns
                    for state_id_2 = 1 : Ns
                        fprintf(file_id, '%0.18e %0.18e\n', real(dissipator(state_id_1, state_id_2)), imag(dissipator(state_id_1, state_id_2)));
                    end
                end
                fclose(file_id);
            end
            
             if(save_super_operator >= 0)
                dissipator_conj = dissipator';
                file_name = sprintf('%sdissipator_conj_%d_%s.txt', data_path, dissipator_id, file_name_suffix);
                file_id = fopen(file_name, 'w');
                for state_id_1 = 1 : Ns
                    for state_id_2 = 1 : Ns
                        fprintf(file_id, '%0.18e %0.18e\n', real(dissipator_conj(state_id_1, state_id_2)), imag(dissipator_conj(state_id_1, state_id_2)));
                    end
                end
                fclose(file_id);
             end
            
             if(save_super_operator >= 0)
                dissipator_mult = dissipator'*dissipator;
                file_name = sprintf('%sdissipator_mult_%d_%s.txt', data_path, dissipator_id, file_name_suffix);
                file_id = fopen(file_name, 'w');
                for state_id_1 = 1 : Ns
                    for state_id_2 = 1 : Ns
                        fprintf(file_id, '%0.18e %0.18e\n', real(dissipator_mult(state_id_1, state_id_2)), imag(dissipator_mult(state_id_1, state_id_2)));
                    end
                end
                fclose(file_id);
            end
            
            super_rp_matrix = super_rp_matrix + ...
                g * 0.5 * (2.0 * kron(eye(Ns),dissipator) * kron(transpose(dissipator'),eye(Ns)) - ...
                kron(transpose(dissipator'*dissipator), eye(Ns)) - ...
                kron(eye(Ns), dissipator'*dissipator));
            
            
            if(save_super_operator >= 0)
                file_name = sprintf('%slindbladian_%s.txt', data_path, file_name_suffix);
                file_id = fopen(file_name, 'w');
                for state_id_1 = 1 : Ns * Ns
                    for state_id_2 = 1 : Ns * Ns
                        fprintf(file_id, '%0.18e %0.18e\n', real(super_rp_matrix(state_id_1, state_id_2)), imag(super_rp_matrix(state_id_1, state_id_2)));
                    end
                end
                fclose(file_id);
            end
        end
        
        if (border_conditions == 1)
            
            dissipator = zeros(Ns);
            for state_id_1 = 1:Ns
                dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), Nc) - bitget(idtox(state_id_1), 1);
                for state_id_2 = 1:Ns
                    if(is_adjacent(idtox(state_id_1), idtox(state_id_2)))
                        hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(state_id_1), idtox(state_id_2)),Nc) == '1' );
                        if(hopping_ids(2) == Nc)
                            if(bitget(idtox(state_id_1), Nc))
                                dissipator(state_id_2,state_id_1) = exp(sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = -exp(-sqrt(-1)*alpha);
                            else
                                dissipator(state_id_2,state_id_1) = -exp(-sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = exp(sqrt(-1)*alpha);
                            end
                        end
                    end
                end
            end
            
            super_rp_matrix = super_rp_matrix + ...
                g * 0.5 * (2.0 * kron(eye(Ns),dissipator) * kron(transpose(dissipator'),eye(Ns)) - ...
                kron(transpose(dissipator'*dissipator), eye(Ns)) - ...
                kron(eye(Ns), dissipator'*dissipator));
        end
        
    elseif(dissipator_type == 0)
        for dissipator_id = 1:Nc
            dissipator=zeros(Ns);
            for state_id_1=1:Ns
                dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), dissipator_id);
            end
            
            super_rp_matrix = super_rp_matrix + ...
                g * 0.5 * (2.0 * kron(eye(Ns),dissipator) * kron(transpose(dissipator'),eye(Ns)) - ...
                kron(transpose(dissipator'*dissipator), eye(Ns)) - ...
                kron(eye(Ns), dissipator'*dissipator));        
        end
    end
    
    if(save_super_operator >= 0)
        file_name = sprintf('%slindbladian_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1 : Ns * Ns
            for state_id_2 = 1 : Ns * Ns
                fprintf(file_id, '%0.18e %0.18e\n', real(super_rp_matrix(state_id_1, state_id_2)), imag(super_rp_matrix(state_id_1, state_id_2)));
            end
        end
        fclose(file_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Zero eigen vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sparse_super_rp_matrix = sparse(super_rp_matrix);
    super_rp_matrix = 0;
    
    [zero_evec, zero_eval] = eigs(sparse_super_rp_matrix, 1, 'sm');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Rho in direct basis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rho = zeros(Ns, Ns);
    for state_id_1 = 1:Ns
        rho(:, state_id_1) = zero_evec(1+(state_id_1-1)*(Ns) : state_id_1*(Ns));
    end
    rho = rho / trace(rho);
    
    if(save_type >= 1)
        file_name = sprintf('%sdiag_rho_in_direct_basis_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            fprintf(file_id, '%0.18e\n', rho(state_id_1, state_id_1));
        end
        fclose(file_id);
    end
    
    if(save_type >= 2)
        file_name = sprintf('%srho_in_direct_basis_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e %0.18e\n', real(rho(state_id_1, state_id_2)), imag(rho(state_id_1, state_id_2)) );
            end
        end
        fclose(file_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Rho in stationary basis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rho_in_stationary_basis = evecs_H' * rho * evecs_H;
    
    if(save_type >= 1)
        file_name = sprintf('%sdiag_rho_in_stationary_basis_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            fprintf(file_id, '%0.18e\n', rho_in_stationary_basis(state_id_1, state_id_1));
        end
        fclose(file_id);
    end
    
    if(save_type >= 2)
        file_name = sprintf('%srho_in_stationary_basis_%s.txt', data_path, file_name_suffix);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e %0.18e\n', real(rho_in_stationary_basis(state_id_1, state_id_2)), imag(rho_in_stationary_basis(state_id_1, state_id_2)) );
            end
        end
        fclose(file_id);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Characteristics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    entropy = -trace(rho * logm(rho)); % entropy is direct basis
    n_part=zeros(1,Nc);
    for state_id_1 = 1:Ns
        n_part = n_part + real(rho(state_id_1, state_id_1)) * (dec2bin(idtox(state_id_1), Nc) == '1');
    end
    imbalance = sum(n_part(1:2:Nc-1) - n_part(2:2:Nc)) / sum(n_part);
    
    file_name = sprintf('%scharacteristics_%s.txt', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    fprintf(file_id, '%0.18e %0.18e', entropy, imbalance);
    fclose(file_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      Checking Discrepancy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    discrepancy = -sqrt(-1) * (H*rho - rho*H);
    
    if(dissipator_type == 1)
        for dissipator_id = 1:Nc-1
            dissipator = zeros(Ns);
            for state_id_1 = 1:Ns
                dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), dissipator_id) - bitget(idtox(state_id_1), dissipator_id+1);
                for state_id_2 = 1:Ns
                    if(is_adjacent(idtox(state_id_1), idtox(state_id_2)))
                        hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(state_id_1), idtox(state_id_2)),Nc) == '1' );
                        if(hopping_ids(2) == dissipator_id)
                            if(bitget(idtox(state_id_1), dissipator_id))
                                dissipator(state_id_2,state_id_1) = exp(sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = -exp(-sqrt(-1)*alpha);
                            else
                                dissipator(state_id_2,state_id_1) = -exp(-sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = exp(sqrt(-1)*alpha);
                            end
                        end
                    end
                end
            end
            
            discrepancy = discrepancy + g * 0.5 * ...
                (2.0 * dissipator * rho * dissipator' - ...
                rho * dissipator' * dissipator - ...
                dissipator' * dissipator * rho);
        end
        
        if (border_conditions == 1)
            
            dissipator = zeros(Ns);
            for state_id_1 = 1:Ns
                dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), Nc) - bitget(idtox(state_id_1), 1);
                for state_id_2 = 1:Ns
                    if(is_adjacent(idtox(state_id_1), idtox(state_id_2)))
                        hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(state_id_1), idtox(state_id_2)),Nc) == '1' );
                        if(hopping_ids(2) == Nc)
                            if(bitget(idtox(state_id_1), Nc))
                                dissipator(state_id_2,state_id_1) = exp(sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = -exp(-sqrt(-1)*alpha);
                            else
                                dissipator(state_id_2,state_id_1) = -exp(-sqrt(-1)*alpha);
                                dissipator(state_id_1,state_id_2) = exp(sqrt(-1)*alpha);
                            end
                        end
                    end
                end
            end
            
            discrepancy = discrepancy + g * 0.5 * ...
                (2.0 * dissipator * rho * dissipator' - ...
                rho * dissipator' * dissipator - ...
                dissipator' * dissipator * rho);
            
        end
    end
    
    if(dissipator_type == 0)
        for dissipator_id = 1:Nc
            dissipator = zeros(Ns);
            for state_id_1 = 1:Ns
                dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), dissipator_id);
            end
            
            discrepancy = discrepancy + g * 0.5 * ...
                (2.0 * dissipator * rho * dissipator' - ...
                rho * dissipator' * dissipator - ...
                dissipator' * dissipator * rho);
        end
    end
    
    error = max(max(abs(discrepancy)));
    
    file_name = sprintf('%schecking_%s.txt', data_path, file_name_suffix);
    file_id = fopen(file_name, 'w');
    fprintf(file_id, '%0.18e %0.18e\n', error, abs(zero_eval));
    fclose(file_id);
    
    toc
    
end