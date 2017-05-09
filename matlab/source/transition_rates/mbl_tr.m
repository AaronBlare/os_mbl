clear all;

filename = 'config.txt';
input_data = importdata(filename);

W = input_data(1); 	% disorder
U = input_data(2);  % interaction
J = input_data(3);  % hopping
g = input_data(4);  % gamma;

Nc = input_data(5); % number of cells

diss_type = input_data(6); 		% Dissipator type: 0-Poletti, 1-Diehl
alpha = input_data(7); 			% Dissipator Diehl phase 
energy_type = input_data(8);   	% 0 if regular, 1 if zero mean

seed = input_data(9);

dump_aux = input_data(10);

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

if(dump_aux >= 1)
    file_name = sprintf('participation_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
    file_id = fopen(file_name, 'w');
    for state_id = 1:Ns
        fprintf(file_id, '%0.18e\n', participation(state_id));
    end
    fclose(file_id);
end

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
        
        transition_rates=transition_rates+abs(Ev'*A*Ev).^2;
    end
    
elseif(diss_type == 0)
    for k=1:Nc
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k);
        end
        
        transition_rates = transition_rates + abs(Ev'*A*Ev).^2;
    end  
end

transition_rates = transition_rates * g;
TR = transition_rates;

for ik=1:Ns
    TR(ik,ik)=TR(ik,ik)-sum(TR(:,ik));
end

if(dump_aux >= 2)
	file_name = sprintf('transition_rates_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
	file_id = fopen(file_name, 'w');
	for state_id_1 = 1:Ns
		for state_id_2 = 1:Ns
			fprintf(file_id, '%0.18e\n', TR(state_id_1, state_id_2));
		end  
	end
	fclose(file_id);
end

[TRev,TReig] = eigs(real(TR),1,'sm');
TRev = TRev / sum(TRev);

if(dump_aux >= 1)
    file_name = sprintf('rho_diag_in_stationary_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
    file_id = fopen(file_name, 'w');
    for state_id = 1:Ns
        fprintf(file_id, '%0.18e\n', TRev(state_id));
    end
    fclose(file_id);
end

rho_in_direct_basis = Ev * diag(TRev) * Ev';

if(dump_aux >= 1)
    file_name = sprintf('rho_diag_in_direct_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
    file_id = fopen(file_name, 'w');
    for state_id = 1:Ns
        fprintf(file_id, '%0.18e\n', rho_in_direct_basis(state_id, state_id));
    end
    fclose(file_id);
end

entropy = -trace(rho_in_direct_basis * logm(rho_in_direct_basis));
n_part=zeros(1,Nc);
for k = 1:Ns
    n_part = n_part + real(rho_in_direct_basis(k,k)) * ( dec2bin(idtox(k),Nc) == '1' );
end
imbalance = sum(n_part(1:2:Nc-1) - n_part(2:2:Nc)) / sum(n_part);

file_name = sprintf('characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
file_id = fopen(file_name, 'w');
fprintf(file_id, '%0.18e %0.18e', entropy, imbalance);
fclose(file_id);

if(dump_aux >= 2)
	file_name = sprintf('rho_in_direct_basis_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', Nc, diss_type, alpha, energy_type, W, U, J, g, seed);
	file_id = fopen(file_name, 'w');
	for state_id_1 = 1:Ns
		for state_id_2 = 1:Ns
			fprintf(file_id, '%0.18e\n', rho_in_direct_basis(state_id_1, state_id_2));
		end  
	end
	fclose(file_id);
end

checking = sum(TRev).^2

toc