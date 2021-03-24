%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
% Number of subsystems
Ni=[5 10 15 25 30 50 55 75 100 125 150 175 200 250 350 500];
alpha0=0.001;

N_h2syn=10; 
N_alg1_sof=7; 
N_alg1_dof=4; 
% subsystem dimensions
nx=3;
nu=1;
ny=1;
C0=[1 0 0];
Vi=0.1;
for i=1:length(Ni)
    N=Ni(i);
    % sparsity patterns
    S=eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
    Sa=kron(S,ones(nx));
    
    A=Sa.*randn(N*nx);
    B=zeros(N*nx,N);
    for j=1:N
        B((j-1)*nx+1:j*nx,j)=[sign(-0.5+rand(1))*(0.1+rand(1)); zeros(nx-1,1)];
    end
    C=kron(eye(N),C0);
    
    A=1.25*A/(max(abs(eig(A)))); 
    Q=eye(nx*N);
    QT=0*Q;
    R=0.1*eye(nu*N);
    
    % generalized plant
    B1=[eye(nx*N) zeros(nx*N,N*ny)];
    C1=[sqrt(Q); zeros(N,nx*N)];
    D11=zeros(size(C1,1), size(B1,2));
    D12=[zeros(nx*N,N); sqrt(R)];
    D21=[zeros(N,nx*N) Vi*eye(N)];
    D22=zeros(N);
    
    Pg=ss(A,[B1 B],[C1 ; C],[D11 D12; D21 D22],-1);
    
    % h2syn:
    if i<=N_h2syn
        tic
        K1=h2syn(Pg,N,N);
        t1(i)=toc
    end
    
    % Algorithm 2:
    Nmax=2000;
    Tmax=50;
    Sp=kron(S,ones(nu,ny));
    
    As=sparse(A);
    Bs=sparse(B);
    Cs=sparse(C);
    Qs=sparse(Q);
    QTs=sparse(QT);
    Rs=sparse(R);
    Ss=sparse(S);
    
    tic 
    K2=alg2_sof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss);
    t2(i)=toc
  
    tic 
    K3=alg2_dof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss,3);
    t3(i)=toc
    
    % Algorithm 1:
    if i<=N_alg1_sof
        tic
        [K4, f] = alg1_sof(Pg,Sp);
        t4(i)=toc
    end
    if i<=N_alg1_dof
        tic
        [K5, f] = alg1_dof(Pg,S,3);
        t5(i)=toc
    end
    
   
end
%%
figure(1)
plot(Ni(1:N_h2syn),t1,'k','LineWidth',1.5)
hold on
grid on
plot(Ni(1:N_alg1_sof),t4,'r','LineWidth',1.5)
plot(Ni,t2,'b','LineWidth',1.5)
plot(Ni(1:N_alg1_dof),t5,'r-.','LineWidth',1.5)
plot(Ni,t3,'b-.','LineWidth',1.5)
legend('h2syn','Alg. 1_S','Alg. 2_S','Alg. 1_D','Alg. 2_D')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
ylim([0 60])