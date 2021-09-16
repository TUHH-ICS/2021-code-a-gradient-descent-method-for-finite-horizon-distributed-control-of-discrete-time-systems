%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all; close all; clc
% Number of subsystems
Ni=[5 10 13 25 30 45 75 100 125 175 250 350 500];
alpha0=0.01;

N_h2syn=9; 
N_systune_sof=6;
N_systune_dof=3; 
N_dof=11;
% subsystem dimensions
nx=3;
nu=1;
ny=1;
C0=[1 0 0];
Vi=0.1;
for l=1:5
for i=1:length(Ni)
    N=Ni(i);
    % sparsity patterns
    flag=1;
    % generate random sparsity patterns
    while(flag)
        St = rand(N) < 5/N; % a logical array consuming little memory
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
    
    Pg=ss(A,[B1 B],[C1 ; C],[D11 D12; D21 D22],1);
    
    % h2syn:
    if i<=N_h2syn
        tic
        K1=h2syn(Pg,N,N);
        t1(l,i)=toc
    end
    
    % Algorithm 1:
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
    
    tic 
    K2=alg1_sof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss);
    t2(l,i)=toc
    
    if i<=N_dof
        tic 
        K3=alg1_dof(As,Bs,Cs,Qs,QTs,Rs,Vi,alpha0,Nmax,Tmax,Ss,3);
        t3(l,i)=toc
    end
    
    % systune:
    if i<=N_systune_sof
        K4=realp('K4',kron(S, zeros(nu,ny)));
        K4.Free=kron(S, ones(nu,ny));
        CL0=lft(Pg,K4);
        CL0.InputName='in';
        CL0.OutputName='out';
        goal=TuningGoal.Variance('in','out',0.01);
        tic
        [CL,fSoft] = systune(CL0,goal);
        t4(l,i)=toc
    end
    if i<=N_systune_dof
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
        tic
        [CL,fSoft] = systune(CL0,goal);
        t5(l,i)=toc
    end  
end
end

%%
t1_=sum(t1)./5;
t2_=sum(t2)./5;
t3_=sum(t3)./5;
t4_=sum(t4)./5;
t5_=sum(t5)./5;
%%
figure(1)
plot(Ni(1:N_h2syn),t1_,'k','LineWidth',1.5)
hold on
grid on
plot(Ni(1:N_systune_sof),t4_,'r','LineWidth',1.5)
plot(Ni,t2_,'b','LineWidth',1.5)
plot(Ni(1:N_systune_dof),t5_,'r-.','LineWidth',1.5)
plot(Ni(1:N_dof),t3_,'b-.','LineWidth',1.5)
legend('h2syn','systune_S','Alg. 1_S','systune_D','Alg. 1_D')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
ylim([0 60])