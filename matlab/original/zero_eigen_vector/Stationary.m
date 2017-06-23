function [A,B] = Stationary

Nc = 8; % size of array
Np = Nc/2; % number of particles
Ns=nchoosek(Nc,Np); % number of states

W=10; % disorder
J=1; % hopping
U=1; % interaction

Diehl=0;
Poletti=1;

tic

num = precalc_states (Nc, Np);
%display ('all values');

% uncomment to view listing of states and adjacency matrix
%
for i = 0:(2 ^ Nc - 1)
    display (sprintf ('x = %2d | bits = %s | is state = %d | id_state = %2d', ...
        i, dec2bin(i, Nc), xtoid(i) > 0, xtoid(i)));
end

for i = 1:num
    display (sprintf ('id = %2d | x = %2d | bits = %s', ...
        i, idtox(i), dec2bin(idtox(i), Nc)));
end
for i = 1:num
    for j = 1:num
        display (sprintf ('(%2d, %s) adj (%2d, %s) = %d', ...
            i,  dec2bin(idtox(i), Nc), j, dec2bin(idtox(j), Nc), ...
            is_adjacent(idtox(i), idtox(j))));
    end
end
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = zeros(Nc); % Halimtonian
Hd=zeros(Nc); % Disorder
Hi=zeros(Nc); % interaction
Hh=zeros(Nc); % hopping

%%%%%%%%%%%%%%
% disorder
%%%%%%%%%%%%%%
rng(1);
E=2*rand(Nc,1)-1;
for k=1:Ns
    Hd(k,k)=(dec2bin(idtox(k), Nc)=='1')*E;
end

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

H=-J*Hh+U*Hi+2*W*Hd;

save 'H.txt' H -ascii
save 'Hd.txt' Hd -ascii
save 'Hi.txt' Hi -ascii
save 'Hh.txt' Hh -ascii

Ns=Ns

[Ev,Eg]=eig(H); % Anderson modes

Ev_original=zeros(Nc,Ns);

for k=1:Ns
    for kk=1:Ns
        Ev_original(:,k)=Ev_original(:,k)+Ev(kk,k)*(dec2bin(idtox(kk), Nc)=='1')';
    end
end

toc

figure(1)
plot(diag(Eg)',sum(Ev.^4),'o-')
xlabel('eig')
ylabel('IPR')

figure(2)
h=pcolor(log10(Ev.^2+1e-8));
set(h,'EdgeColor','None');
colorbar
xlabel('eigenvector No')
ylabel('eigenvector element')


figure(3)
h=pcolor((Ev_original.^1+0e-8));
set(h,'EdgeColor','None');
colorbar
xlabel('eigenvector number')
ylabel('element in occupation number space')

figure(4)
plot(diag(Eg)',sum(diff(Ev_original).^2),'o-')
xlabel('eig')
ylabel('next neighbor difference in occupation')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(Diehl)
    
    for k=1:Nc-1
        
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k)-bitget(idtox(kk),k+1);
            for kkk=1:Ns
                if(is_adjacent(idtox(kk), idtox(kkk)))
                    hop=1+Nc-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),Nc)=='1');
                    if(hop(2)==k)
                        
                        if((bitget(idtox(kk),k)))
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
        
        save 'A.txt' A -ascii -append
        
    end
    
end


if(Poletti)
    
    for k=1:Nc
        
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k);
        end
        
        
        save 'A.txt' A -ascii -append
        
    end
end

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
