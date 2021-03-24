%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
% Number of subsystems
Ni=[5 10 15 30 50 75 100 125 150 175 200 250 350 500];
alpha0=0.001;
N_alg1=5; 
N_lmi=3;  
N_dlqr=12; 

% subsystem dimensions
nx=3;
nu=1;
%%
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
    A=1.5*A/max(abs(eig(A))); 
    Q=eye(nx*N);
    QT=0*Q;
    R=0.1*eye(nu*N);
    if i<=N_dlqr
        tic
        K1=dlqr(A,B,Q,R);
        t1(i)=toc
    end
    
    Sp=kron(S,ones(nu,nx));
    Tmax=25; 
    Nmax=1000;
    As=sparse(A);
    Bs=sparse(B);
    Qs=sparse(Q);
    QTs=sparse(QT);
    Rs=sparse(R);
    Ss=sparse(Sp);
    tic 
    K2=alg2_sf(As,Bs,Qs,QTs,Rs,Nmax,Tmax,Ss,alpha0);
    t2(i)=toc
    
    if i<=N_alg1
        Qx=Q;
        Qu=R;
        Qxu=zeros(N*nx,N*nu);
        tic
        K3=alg1_sf(A,B,Qx,Qxu,Qu,Sp);
        t3(i)=toc
    end
    
    if i<=N_lmi
        P.A=A;
        P.B1=eye(nx*N);
        P.B2=B;
        P.C1=[sqrt(Q); zeros(N,nx*N)];
        P.D12=[zeros(nx*N,N); sqrt(R)];
        tic
        K4=LMI_sf(P,N,Sp);
        t4(i)=toc
    end

end

figure(1)
plot(Ni(1:N_dlqr),t1,'k','LineWidth',1.5)
hold on
grid on
plot(Ni(1:N_alg1),t3,'r','LineWidth',1.5)
plot(Ni,t2,'b','LineWidth',1.5)
plot(Ni(1:N_lmi),t4,'color',[0 0.5 0],'LineWidth',1.5)
legend('dlqr','Alg. 1','Alg. 2','LMI')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
ylim([0 60])
