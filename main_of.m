%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
NoS=50;
N=10;
C0=[1 0 0];
C=kron(eye(N),C0);
%dimensions of subsystems
nx=3;
nu=1;
ny=1;
% sparsity patterns
S=eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
Sa=kron(S,ones(nx));

Vi=0.1;
for i=1:NoS
    A=Sa.*randn(N*nx);
    B=zeros(N*nx,N);
    for j=1:N
        B((j-1)*nx+1:j*nx,j)=[sign(-0.5+rand(1))*(0.1+rand(1)); zeros(nx-1,1)];
    end
    
    A=1.25*A/(max(abs(eig(A)))); 
    Q=eye(nx*N);
    QT=0*Q;
    R=0.1*eye(nu*N);
    
    % generalized plant
    B1=[eye(nx*N) zeros(nx*N,N*ny)];
    C1=[sqrt(Q); zeros(N,nx*N)];
    D11=zeros(size(C1,1), size(B1,2));
    D12=[zeros(nx*N,N); sqrt(R)];
    D21=[zeros(N,nx*N) sqrt(Vi)*eye(N)];
    D22=zeros(N);
    
    Pg=ss(A,[B1 B],[C1 ; C],[D11 D12; D21 D22],-1);
    
    % h2syn:
    K1=h2syn(Pg,N,N);
    
    % Algorithm 2:
    alpha0=0.001;
    Nmax=2000; 
    Tmax=50;    
    Sp=kron(S,ones(nu,ny));
    
    As=sparse(A);
    Bs=sparse(B);
    Cs=sparse(C);
    Qs=sparse(Q);
    QTs=sparse(QT);
    Rs=sparse(R);
    Ss=sparse(Sp);
    
    % static output feedback
    K2=alg2_sof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss);
    
    % dynmaic output feedback
    nkx=3;
    K3=alg2_dof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss,nkx);

    % Algorithm 1:
    % static output feedback
    K4=alg1_sof(Pg,Sp);
    
    % dynamic output feedback
    nkx=3;
    K5=alg1_dof(Pg,S,nkx);
    
    N1(i)=norm(lft(Pg,K1));

    N2(i)=norm(lft(Pg,K2));
    
    N3(i)=norm(lft(Pg,K3));
    
    N4(i)=norm(lft(Pg,K4));
    
    N5(i)=norm(lft(Pg,K5));
   
end
%%
P1=N1./N2;
P2=N1./N3;
P3=N1./N4;
P4=N1./N5;

P1_=fliplr(sort([P1 1]));
P2_=fliplr(sort([P2 1]));
P3_=fliplr(sort([P3 1]));
P4_=fliplr(sort([P4 1]));

x_=0:(1/NoS):1;

figure(1)
plot(100*x_,P3_,'r','LineWidth',1.5)
hold on
grid on
box on
plot(100*x_,P1_,'b','LineWidth',1.5)
plot(100*x_,P4_,'r-.','LineWidth',1.5)
plot(100*x_,P2_,'b-.','LineWidth',1.5)
legend('Alg. 1_S','Alg. 2_S','Alg. 1_D','Alg. 2_D')
xlabel('Fraction of systems [%]')

ylabel('Centralized/ distribtuted')
ylim([0 1])
