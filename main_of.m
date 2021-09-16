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
   
    B=zeros(N*nx,N);
    
    flag=1;
    % generate random sparsity patterns
    while(flag)
        St = rand(N) < 0.3; % a logical array consuming little memory
        S = (tril(St)+tril(St)'+eye(N))>0;
        G=graph(S);
        bins = conncomp(G);
        % check if graph is connected
        if bins==ones(1,N)
            flag=0;
        else
            flag=1;
        end
    end
    Sa=kron(S,ones(nx));
    
    % number of neighbors
    count_n=-ones(1,N);
    for j=1:N
        for k=1:N
            if S(j,k)==1
                count_n(j)=count_n(j)+1;
            end
        end
    end
    
    A=Sa.*randn(N*nx);
    
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
    
    Pg=ss(A,[B1 B],[C1 ; C],[D11 D12; D21 D22],1);
    
    % h2syn:
    K1=h2syn(Pg,N,N);
    
    % Algorithm 1:
    alpha0=0.01;
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
    K2=alg1_sof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss);
    
    % dynmaic output feedback
    nkx=3;
    K3=alg1_dof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss,nkx);

    % systune:
    % static output feedback
    K4=realp('K4',kron(S, zeros(nu,ny)));
    K4.Free=kron(S, ones(nu,ny));
    CL0=lft(Pg,K4);
    CL0.InputName='in';
    CL0.OutputName='out';
    goal=TuningGoal.Variance('in','out',0.01);
    [CL,fSoft] = systune(CL0,goal);
    K4=CL.A.Blocks.K4.value;
    
    % dynamic output feedback
    nkx=3;
    AK=realp('AK',kron(S, 0.1*randn(nkx,nkx)));
    AK.Free=kron(S, ones(nkx,nkx));
    BK=realp('BK',kron(eye(N), 0.1*randn(nkx,ny)));
    BK.Free=kron(eye(N), ones(nkx,ny));
    CK=realp('CK',kron(S, 0.1*randn(nu,nkx)));
    CK.Free=kron(S, ones(nu,nkx));
    DK=realp('DK',kron(eye(N), 0.1*randn(nu,ny)));
    DK.Free=kron(eye(N), ones(nu,ny));
    K5=ss(AK,BK,CK,DK,1);
    CL0=lft(Pg,K5);
    CL0.InputName='in';
    CL0.OutputName='out';
    goal=TuningGoal.Variance('in','out',0.01);
    [CL,fSoft] = systune(CL0,goal);
    AK=CL.A.Blocks.AK.value;
    BK=CL.A.Blocks.BK.value;
    CK=CL.A.Blocks.CK.value;
    DK=CL.A.Blocks.DK.value;
    K5=ss(AK,BK,CK,DK,1);

    
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
legend('systune_S','Alg. 1_S','systune_D','Alg. 1_D')
xlabel('Fraction of systems [%]')

ylabel('Centralized/ distribtuted')
ylim([0 1])
