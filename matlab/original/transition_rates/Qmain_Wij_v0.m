%Markovian calculations

function Qmain_SVD

global N

show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);
periodic=0; % periodic b.c.
dissipator_kind=1; % put 1 if symmetric, -1 if antisymmetric

seed = 1;

Diehl=1;
Poletti=0;

tic

N = 10; % size of array
NN = N/2; % number of particles
Ns=nchoosek(N,NN); % number of states

J=1; % hopping constant
W=8; % disorder
U=2; % on-site interaction
g=0.1; % gamma;

Hopping=zeros(Ns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = precalc_states (N, NN);


H = zeros(N); % Halimtonian
Hd=zeros(N); % Disorder
Hi=zeros(N); % interaction
Hh=zeros(N); % hopping

%%%%%%%%%%%%%%
% disorder
%%%%%%%%%%%%%%

rng(seed);

E=2*rand(N,1)-1;
%E=E-sum(E)/N;

B=zeros(N,Ns);

for k=1:Ns
    Hd(k,k)=(1-(dec2bin(idtox(k), N)=='1'))*E;
    B(:,k)=(1-(dec2bin(idtox(k), N)=='1'))';
end

save 'B.txt' B -ascii

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

H=-J*Hh+U*Hi+W*Hd;

save 'H.txt' H -ascii
save 'Hd.txt' Hd -ascii
save 'Hi.txt' Hi -ascii
save 'Hh.txt' Hh -ascii

Ns=Ns

[Ev,Eg]=eig(H); % Anderson modes

Participation=1./sum(abs(Ev).^4,1);
neighbor_hopping=diag(Ev'*(-Hh)*Ev);

%B=B
%Occupation=B*abs(Ev).^2;

figure;
subplot(1,2,1)
plot(Participation,'o-')
xlabel('q')
ylabel('Participation')
subplot(1,2,2)
plot(neighbor_hopping,'o-')
xlabel('q')
ylabel('Hh')

%figure;
%hold on
%plot(Occupation(:,1),'o-')
%plot(Occupation(:,2),'o-')
%plot(Occupation(:,3),'o-')
%plot(Occupation(:,4),'o-')
%plot(Occupation(:,5),'o-')

%figure;
%plot(Ev(:,1),'o-')
%hold on
%plot(Ev(:,2),'o-')
%plot(Ev(:,3),'o-')

%figure;
%hold on
%plot(Occupation(:,Ns-2),'o-')
%plot(Occupation(:,Ns-1),'o-')
%plot(Occupation(:,Ns),'o-')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(Diehl)
    
    for k=1:N-1
        
        
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k)-bitget(idtox(kk),k+1);
            for kkk=1:Ns
                if(is_adjacent(idtox(kk), idtox(kkk)))
                    hop=1+N-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),N)=='1');
                    if(hop(2)==k)
                        
                        if((bitget(idtox(kk),k)))
                            A(kkk,kk)=1*dissipator_kind;
                            A(kk,kkk)=-1*dissipator_kind;
                        else
                            A(kkk,kk)=-1*dissipator_kind;
                            A(kk,kkk)=1*dissipator_kind;
                        end
                        
                    end
                end
                
            end
        end
        
        Hopping=Hopping+abs(Ev'*A*Ev).^2;
    end
    
end


if(Poletti)
    
    for k=1:N
        
        
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k);
        end
        
        Hopping=Hopping+abs(Ev'*A*Ev).^2;
        
        
    end
    
end

Hopping=Hopping*g;
M=Hopping;


for ik=1:Ns
    M(ik,ik)=M(ik,ik)-sum(M(:,ik));
end
[Mev,Meig]=eigs(real(M),1,'sm');
Mev=Mev/sum(Mev);

%save 'M.txt' M -ascii


Rho=Ev*diag(Mev)*Ev';


figure;
plot(diag(M),'o-')
xlabel('q')
ylabel('M(q,q)')


MA=M-diag(diag(M));
%M=expm(4*M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc



Hop_plot=padarray(MA,[1 1],'post');
%  M_plot=padarray(M,[1 1],'post');

Rho_plot=padarray(Rho,[1 1],'post');


figure;
[nx,ny]=meshgrid([0:Ns],[0:Ns]);
h=pcolor(nx,ny,(Hop_plot));
set(h,'EdgeColor','None')
xlabel('q')
ylabel('p')
title('Hop_{q,p}')
colormap('hot')
colorbar


% figure;
% [nx,ny]=meshgrid([0:Ns],[0:Ns]);
% h=pcolor(nx,ny,M_plot);
% set(h,'EdgeColor','None')
% xlabel('q')
% ylabel('p')
% title('M_{q,p}')
% colormap('hot')
% colorbar

%modularity_matrix_plot(MA);


figure;
hold on
plot(Mev,'s-r')
%plot(mass_center,diag(abs(Rho_and)),'^g')
xlabel('q')
ylabel('|\rho_{q,q}|')

figure;
plot(neighbor_hopping,Mev,'or')

xlabel('Hh')
ylabel('|\rho_{q,q}|')

figure;
[nx,ny]=meshgrid([0:Ns],[0:Ns]);
h=pcolor(nx,ny,abs(Rho_plot));
set(h,'EdgeColor','None')
xlabel('m')
ylabel('n')
title('abs(\rho_{m,n})')
colormap('hot')
colorbar


figure;
plot(diag(Eg),'o-')
xlabel('q')
ylabel('\epsilon(q)')

InvPart_and=sum(Mev).^2

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
