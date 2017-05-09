clear all;

filename = 'config_many.txt';
input_data = importdata(filename);

W = input_data(1); 	% disorder
U = input_data(2);  % interaction
J = input_data(3);  % hopping
g = input_data(4);  % gamma;

Nc = input_data(5); % number of cells

diss_type = input_data(6); 		% Dissipator type: 0-Poletti, 1-Diehl
alpha = input_data(7); 			% Dissipator Diehl phase
energy_type = input_data(8);   	% 0 if regular, 1 if zero mean

seed_start = input_data(9);
seed_num = input_data(10);

dump_aux = input_data(11);

for seed = seed_start : seed_start + (seed_num - 1)
    
    Np = Nc/2; % number of particles
    Ns = nchoosek(Nc,Np); % number of states
    init = precalc_states(Nc, Np);
    
    tic
    
    transition_rates = zeros(Ns);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating H
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    H  = zeros(Ns); % Halimtonian
    Hd = zeros(Ns); % Disorder
    Hi = zeros(Ns); % Interaction
    Hh = zeros(Ns); % Hopping
    
    %%%%%%%%%%%%%%
    % Disorder
    %%%%%%%%%%%%%%
    rng(seed)
    E = 2 * rand(Nc,1) - 1;
    
    if energy_type == 1
        E = E - sum(E) / Nc;
    end
    
    
    for k = 1:Ns
        Hd(k,k) = ( dec2bin(idtox(k), Nc) == '1' ) * E;
    end
    
    file_name = sprintf('random_energies_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
    file_id = fopen(file_name, 'w');
    for dump_id = 1:Nc
        fprintf(file_id, '%0.18e\n', E(dump_id));
    end
    fclose(file_id);
    
    %%%%%%%%%%%%%%%%%%%%%
    % Interaction
    %%%%%%%%%%%%%%%%%%%%%
    for k=1:Ns
        Hi(k,k)=sum(dec2bin(bitand(idtox(k), bitshift(idtox(k), 1))) == '1');
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    % Hopping
    %%%%%%%%%%%%%%%%%%%%%
    
    for k=1:Ns
        for kk=1:Ns
            Hh(k,kk)=is_adjacent(idtox(k), idtox(kk));
        end
    end
    
    %%%%%%%%%%%%
    % Total
    %%%%%%%%%%%%
    
    H = -J*Hh + U*Hi + 2.0*W*Hd;
    
    [Ev,Eg] = eig(H); % Anderson modes
    
    Eg_diag = diag(Eg);
    
    if(dump_aux >= 1)
        file_name = sprintf('eg_hamiltonian_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id = 1:Ns
            fprintf(file_id, '%0.18e\n', Eg_diag(state_id));
        end
        fclose(file_id);
    end
    
    if(dump_aux >= 2)
        file_name = sprintf('ev_hamiltonian_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e\n', Ev(state_id_1, state_id_2));
            end
        end
        fclose(file_id);
    end
    
    participation = 1./sum(abs(Ev).^4,1);
    
    if(dump_aux >= 2)
        file_name = sprintf('participation_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id = 1:Ns
            fprintf(file_id, '%0.18e\n', participation(state_id));
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
    
    if(diss_type == 1)
        
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
        
    elseif(diss_type == 0)
        for k=1:Nc
            A=zeros(Ns);
            for kk=1:Ns
                A(kk,kk)=bitget(idtox(kk),k);
            end
            
            P=P+...
                g/2*(2*kron(eye(Ns),A)*kron(transpose(A'),eye(Ns))-...
                kron(transpose(A'*A),eye(Ns))-kron(eye(Ns),A'*A));
        end
    end
    
    Ps = sparse(P);
    P=0;
    
    [zero_evec,zero_eval]=eigs(Ps,1,'sm');
    stationary_rho_array = zero_evec;
    
    st_rho = zeros(Ns, Ns);
    for i=1:Ns
        st_rho(:,i)=stationary_rho_array(1+(i-1)*(Ns):i*(Ns));
    end
    st_rho=st_rho/trace(st_rho);
    
    if(dump_aux >= 1)
        file_name = sprintf('rho_diag_in_direct_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id = 1:Ns
            fprintf(file_id, '%0.18e\n', st_rho(state_id, state_id));
        end
        fclose(file_id);
    end
    
    if(dump_aux >= 2)
        file_name = sprintf('rho_in_direct_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e %0.18e\n', real(st_rho(state_id_1, state_id_2)), imag(st_rho(state_id_1, state_id_2)) );
            end
        end
        fclose(file_id);
    end
    
    st_rho_and=Ev'*st_rho*Ev;
    
    if(dump_aux >= 1)
        file_name = sprintf('rho_diag_in_stationary_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id = 1:Ns
            fprintf(file_id, '%0.18e\n', st_rho_and(state_id, state_id));
        end
        fclose(file_id);
    end
    
    if(dump_aux >= 2)
        file_name = sprintf('rho_in_stationary_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
        file_id = fopen(file_name, 'w');
        for state_id_1 = 1:Ns
            for state_id_2 = 1:Ns
                fprintf(file_id, '%0.18e %0.18e\n', real(st_rho_and(state_id_1, state_id_2)), imag(st_rho_and(state_id_1, state_id_2)) );
            end
        end
        fclose(file_id);
    end
    
    stationary_entropy = -trace(st_rho * logm(st_rho)); % entropy is direct basis
    n_part=zeros(1,Nc);
    for k=1:Ns
        n_part=n_part+real(st_rho(k,k))*(dec2bin(idtox(k),Nc)=='1');
    end
    stationary_imbalance = sum(n_part(1:2:Nc-1)-n_part(2:2:Nc))/sum(n_part);
    
    file_name = sprintf('characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
    file_id = fopen(file_name, 'w');
    fprintf(file_id, '%0.18e %0.18e', stationary_entropy, stationary_imbalance);
    fclose(file_id);
    
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
    if(diss_type == 1)
        for k = 1:Nc-1
            A = zeros(Ns);
            for kk = 1:Ns
                A(kk,kk) = bitget(idtox(kk),k) - bitget(idtox(kk),k+1);
                for kkk = 1:Ns
                    if(is_adjacent(idtox(kk), idtox(kkk)))
                        hop=1+Nc-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),Nc)=='1');
                        if(hop(2)==k)
                            
                            if((~bitget(idtox(kk),k))&&(bitget(idtox(kkk),k+1)))
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
            
            dy=dy+g/2*(2*A*st_rho*A'-st_rho*A'*A-A'*A*st_rho);
        end
    end
    
    if(diss_type == 0)
        for k = 1:Nc
            A = zeros(Ns);
            for kk = 1:Ns
                A(kk,kk) = bitget(idtox(kk),k);
            end
            
            dy=dy+g/2*(2*A*st_rho*A'-st_rho*A'*A-A'*A*st_rho);
        end
    end
    
    err = max(max(abs(dy)));
    
    file_name = sprintf('info_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
    file_id = fopen(file_name, 'w');
    fprintf(file_id, '%0.18e %0.18e\n', err, abs(zero_eval));
    fclose(file_id);
    
    toc
    
end