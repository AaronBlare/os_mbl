clear all;


for seed = 1:100
    
    seed = seed
    
    image1 = sqrt(-1);
    
    period = 0.1;
    num_period_segments = 1;
    deep = 16;
    
    Nc = 8;                 % number of cells
    Np = Nc/2;              % number of particles
    Ns = nchoosek(Nc, Np);  % number of states
    
    
    W = 10.0; % disorder
    U = 1.0;  % interaction
    J = 1.0;  % hopping
    g = 0.1;  % gamma;
    
    diss_type = 1; % 1-Poletti  2-Diehl
    
    if diss_type == 1
        dis_count = Nc;
    elseif diss_type == 2
        dis_count = Nc - 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating H
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    num = precalc_states(Nc, Np);
    
    H = zeros(Ns);  % Halimtonian
    Hd = zeros(Ns); % Disorder
    Hi = zeros(Ns); % interaction
    Hh = zeros(Ns); % hopping
    
    %%%%%%%%%%%%%%
    % disorder
    %%%%%%%%%%%%%%
    rng(seed)
    E = 2 * rand(Nc,1) - 1;
    for k = 1:Ns
        Hd(k,k) = ( dec2bin(idtox(k), Nc) == '1' ) * E;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    % Interaction
    %%%%%%%%%%%%%%%%%%%%%
    for k = 1:Ns
        Hi(k,k) = sum(dec2bin( bitand(idtox(k), bitshift(idtox(k), 1)) ) == '1');
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    % Hopping
    %%%%%%%%%%%%%%%%%%%%%
    for k = 1:Ns
        for kk = 1:Ns
            Hh(k,kk) = is_adjacent(idtox(k), idtox(kk));
        end
    end
    
    %%%%%%%%%%%%
    % Total H
    %%%%%%%%%%%%
    H = -J*Hh + U*Hi + 2*W*Hd;
    
    [Ev,Eg] = eig(H);
    Ev_t = Ev';
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % aux data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file_name = sprintf('aux_data_period%0.2f_dt%d_Ns%d_W%0.2f_U%0.2f_J%0.2f_g%0.2f_seed%d.bin', period, diss_type, Ns, W, U, J, g, seed);
    file_id = fopen(file_name, 'wb');
    for state_id = 1:Ns
        state_in_bin = idtox(state_id);
        fwrite(file_id, state_in_bin, 'int');
    end
    
    res_H = zeros(2*Ns*Ns,1);
    cur_id = 1;
    for state_id_1 = 1:Ns
        for state_id_2 = 1:Ns
            res_H(cur_id)= real(H(state_id_1, state_id_2));
            res_H(cur_id+1)= imag(H(state_id_1, state_id_2));
            cur_id = cur_id + 2;
        end
    end
    fwrite(file_id, res_H, 'double');
    
    res_Ev = zeros(2*Ns*Ns,1);
    cur_id = 1;
    for state_id_1 = 1:Ns
        for state_id_2 = 1:Ns
            res_Ev(cur_id)= real(Ev(state_id_1, state_id_2));
            res_Ev(cur_id+1)= imag(Ev(state_id_1, state_id_2));
            cur_id = cur_id + 2;
        end
    end
    fwrite(file_id, res_Ev, 'double');
    
    
    res_Ev_t = zeros(2*Ns*Ns,1);
    cur_id = 1;
    for state_id_1 = 1:Ns
        for state_id_2 = 1:Ns
            res_Ev_t(cur_id)= real(Ev_t(state_id_1, state_id_2));
            res_Ev_t(cur_id+1)= imag(Ev_t(state_id_1, state_id_2));
            cur_id = cur_id + 2;
        end
    end
    fwrite(file_id, res_Ev_t, 'double');
    
    fclose(file_id);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dissipator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if(diss_type == 1)      % Poletti
        for k=1:dis_count
            A{k}=zeros(Ns);
            for kk=1:Ns
                A{k}(kk,kk)=bitget(idtox(kk),k);
            end
        end
    elseif(diss_type == 2)  % Diehl
        for k=1:dis_count
            A{k}=zeros(Ns);
            for kk=1:Ns
                A{k}(kk,kk)=bitget(idtox(kk),k)-bitget(idtox(kk),k+1);
                for kkk=1:Ns
                    if(is_adjacent(idtox(kk), idtox(kkk)))
                        hop=1+Nc-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),Nc)=='1');
                        if(hop(2)==k)
                            
                            if((bitget(idtox(kk),k)))
                                A{k}(kkk,kk)=1;
                                A{k}(kk,kkk)=-1;
                            else
                                A{k}(kkk,kk)=-1;
                                A{k}(kk,kkk)=1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    H1=H;
    for l=1:dis_count
        H1 = H1 - image1/2 * (g *((A{l})') * A{l});
    end
    
    res = zeros(2*Ns*Ns,1);
    
    file_name_data = sprintf('input_data_period%0.2f_dt%d_Ns%d_W%0.2f_U%0.2f_J%0.2f_g%0.2f_seed%d.bin', period, diss_type, Ns, W, U, J, g, seed);
    fid = fopen(file_name_data, 'wb');
    
    fwrite(fid, period, 'double');
    fwrite(fid, Ns, 'int');
    fwrite(fid, num_period_segments, 'int');
    fwrite(fid, deep, 'int');
    
    dt = period;
    split = zeros(deep, 1);
    split(:) = 2;
    
    for j = 1:deep
        
        dt = dt/split(j);
        
        fwrite(fid, dt, 'double');
        fwrite(fid, split(j), 'int');
        
        G1 = expm(-image1*(H1)*dt);
        
        m=1;
        for l = 1:Ns
            for k = 1:Ns
                res(m) = real(G1(l,k));
                res(m+1) = imag(G1(l,k));
                m = m + 2;
            end
        end
        
        fwrite(fid, res, 'double');
    end
    
    fwrite(fid, dis_count, 'int');
    
    for j = 1:dis_count
        
        fwrite(fid, g, 'double');
        
        m=1;
        for l = 1:Ns
            for k = 1:Ns
                res(m)= real(A{j}(l,k));
                res(m+1)= imag(A{j}(l,k));
                m = m + 2;
            end
        end
        fwrite(fid, res, 'double');
    end
    
    fclose(fid);
    
end

% create enumeration of n bit states with m ones
function num = precalc_states(n, m)
global states_id
k = 1;
for i = 0:(2 ^ n - 1)
    if ((sum(dec2bin(i) == '1') == 2) && (sum(dec2bin(bitand(i, bitshift(i, 1))) == '1') == 1))
        states_id.adjacent(i + 1) = 1;
    else
        states_id.adjacent(i + 1) = 0;
    end
    if (sum(dec2bin(i) == '1') == m)
        states_id.xtoid(i + 1) = k;
        states_id.idtox(k) = i;
        k = k + 1;
    else
        states_id.xtoid(i + 1) = 0;
    end
end
num = k - 1;
end

% x - allowable bit states (int)
function [id] = xtoid (x)
global states_id
id = states_id.xtoid(x + 1);
end

% id - allowable id (int)
function [x] = idtox (id)
global states_id
x = states_id.idtox(id);
end

% x, y - allowable bit states (int)
% if x is adjacent to y then adj = 1 else adj = 0
function [adj] = is_adjacent (x, y)
global states_id
adj = states_id.adjacent(bitxor(x, y) + 1);
end

