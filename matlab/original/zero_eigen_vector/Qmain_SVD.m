function Qmain_SVD

global N

show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);
periodic=0; % periodic b.c.

Diehl=1;
Poletti=0;

tic

N = 8; % size of array
NN = N/2; % number of particles
Ns=nchoosek(N,NN); % number of states

J=1; % hopping constant
W=10; % disorder
U=1; % on-site interaction
g=0.1; % gamma;

seed = 1;

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
for k=1:Ns
    Hd(k,k)=(dec2bin(idtox(k), N)=='1')*E;
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

%%%%%%%%%%%%%%%%%%%
% Creating matrix P*Rho=0
%%%%%%%%%%%%%%%%%%%
%-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho)=0

P=zeros((Ns)^2);
P=-sqrt(-1)*(kron(eye(Ns),H)-kron(transpose(H),eye(Ns)));


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
        
        
        P=P+...
            g/2*(2*kron(eye(Ns),A)*kron(transpose(A'),eye(Ns))-...
            kron(transpose(A'*A),eye(Ns))-kron(eye(Ns),A'*A));
        
    end
    
end


if(Poletti)
    
    for k=1:N
        
        
        A=zeros(Ns);
        for kk=1:Ns
            A(kk,kk)=bitget(idtox(kk),k);
        end
        
        
        P=P+...
            g/2*(2*kron(eye(Ns),A)*kron(transpose(A'),eye(Ns))-...
            kron(transpose(A'*A),eye(Ns))-kron(eye(Ns),A'*A));
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps=sparse(P);
P=0;
[BReig,Seig]=eigs(Ps,1,'sm');

%[BL,S,BR]=svd(P);
%ii=0;
%for i=1:(N)^2
%    if(abs(S(i,i))<1e-8)
%        ii=i;
%    end
%end

Rho2eig=BReig;
Seig=Seig
%Rho2=BR(:,ii);

%rhoerr=max(max(abs(Rho2eig-Rho2)))

for i=1:Ns
    Rho(:,i)=Rho2eig(1+(i-1)*(Ns):i*(Ns));
end

Rho=Rho/trace(Rho);

toc

%dy=-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho);
%err=max(max(abs(dy)))


Rho_and=Ev'*Rho*Ev;

dRho_and=diag(Rho_and);
[sorted,Rindex]=sort(dRho_and,'descend');

Rho_plot=padarray(Rho,[1 1],'post');

Rho_and_plot=padarray(Rho_and,[1 1],'post');


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
%[nx,ny]=meshgrid([0:N],[0:N]);
h=pcolor(nx,ny,angle(Rho_plot));
set(h,'EdgeColor','None')
xlabel('m')
ylabel('n')
title('arg(\rho_{m,n})')
colormap('hot')
colorbar


figure;
[nx,ny]=meshgrid([0:Ns],[0:Ns]);
h=pcolor(nx,ny,abs(Rho_and_plot));
set(h,'EdgeColor','None')
xlabel('q')
ylabel('p')
title('abs(\rho_{q,p})')
colormap('hot')
colorbar


figure;
%[nx,ny]=meshgrid([0:N],[0:N]);
h=pcolor(nx,ny,angle(Rho_and_plot));
set(h,'EdgeColor','None')
xlabel('q')
ylabel('p')
title('arg(\rho_{q,p})')
colormap('hot')
colorbar

figure;
plot(diag(Eg),diag(abs(Rho_and)),'o-')
xlabel('eig')
ylabel('|\rho_{q,q}|')


figure;
plot(Rindex,'o-')

figure;
subplot(4,1,1)
plot(Ev(:,Rindex(1)),'o-b')
hold on
plot(Ev(:,Rindex(2)),'s-r')
plot(Ev(:,Rindex(3)),'^-g')

subplot(4,1,2)
[mval mindex]=max(abs(Ev(:,Rindex(1))))
plot(dec2bin(idtox(mindex),N)=='1','o-b')
hold on
subplot(4,1,3)
[mval mindex]=max(abs(Ev(:,Rindex(2))));
plot(dec2bin(idtox(mindex),N)=='1','s-r')
subplot(4,1,4)
[mval mindex]=max(abs(Ev(:,Rindex(3))));
plot(dec2bin(idtox(mindex),N)=='1','^-g')

figure;
subplot(4,1,1)
plot(Ev(:,Rindex(Ns)),'o-b')
hold on
plot(Ev(:,Rindex(Ns-1)),'s-r')
plot(Ev(:,Rindex(Ns-2)),'^-g')

subplot(4,1,2)
[mval mindex]=max(abs(Ev(:,Rindex(Ns))))
plot(dec2bin(idtox(mindex),N)=='1','o-b')
hold on
subplot(4,1,3)
[mval mindex]=max(abs(Ev(:,Rindex(Ns-1))));
plot(dec2bin(idtox(mindex),N)=='1','s-r')
subplot(4,1,4)
[mval mindex]=max(abs(Ev(:,Rindex(Ns-2))));
plot(dec2bin(idtox(mindex),N)=='1','^-g')


figure;
plot(diag(Eg),'o-')
xlabel('q')
ylabel('\epsilon(q)')



S=-trace(Rho*logm(Rho)) % entropy is direct basis
S_and=-trace(Rho_and*logm(Rho_and)) % entropy in stationary basis

Npart=zeros(1,N);
for k=1:Ns
    Npart=Npart+real(Rho(k,k))*(dec2bin(idtox(k),N)=='1');
end

Imb=sum(Npart(1:2:N-1)-Npart(2:2:N))/sum(Npart)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=-sqrt(-1)*(H*Rho-Rho*H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Dissipators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:N-1
    
    
    A=zeros(Ns);
    for kk=1:Ns
        A(kk,kk)=bitget(idtox(kk),k)-bitget(idtox(kk),k+1);
        for kkk=1:Ns
            if(is_adjacent(idtox(kk), idtox(kkk)))
                hop=1+N-find(dec2bin(bitxor(idtox(kk),idtox(kkk)),N)=='1');
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
    
    dy=dy+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho);
    
end

err=max(max(abs(dy)))

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
