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
N_systune=5; 
N_lmi=3; 
N_dlqr=12;

% subsystem dimensions
nx=3;
nu=1;
%%
for i=1:length(Ni)
    N=Ni(i);
    % sparsity patterns
    flag=1;
    % generate random sparsity patterns
    while(flag)
        St = rand(N) < 5/N; 
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
    A=1.5*A/max(abs(eig(A))); 
    Q=eye(nx*N);
    QT=0*Q;
    R=0.1*eye(nu*N);
    if i<=N_dlqr
        tic
        K1=dlqr(A,B,Q,R);
        t1(i)=toc
    end
    
    % generalized plant
    Ts=1;
    gen_plant=ss(A,[eye(nx*N) B],[sqrt(Q); zeros(N,nx*N); eye(nx*N)],[zeros(nx*N) zeros(nx*N,N); zeros(N,nx*N) sqrt(R); zeros(nx*N,(nx+1)*N)],Ts);
    
    Sp=kron(S,ones(nu,nx));
    Tmax=25; 
    Nmax=1500;
    As=sparse(A);
    Bs=sparse(B);
    Qs=sparse(Q);
    QTs=sparse(QT);
    Rs=sparse(R);
    Ss=sparse(Sp);
    tic 
    K2=alg1_sf(As,Bs,Qs,QTs,Rs,Nmax,Tmax,Ss,alpha0);
    t2(i)=toc
    
    if i<=N_systune
        Qx=Q;
        Qu=R;
        Qxu=zeros(N*nx,N*nu);
        tic
        K3=realp('K3',kron(S, ones(nu,nx)));
        K3.Free=kron(S, ones(nu,nx));
        CL0=lft(gen_plant,K3);
        CL0.InputName='in';
        CL0.OutputName='out';
        goal=TuningGoal.Variance('in','out',0.01);
        [CL,fSoft] = systune(CL0,goal);
        K3=CL.A.Blocks.K3.value;
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
%%
figure(1)
plot(Ni(1:N_dlqr),t1,'k','LineWidth',1.5)
hold on
grid on
plot(Ni(1:N_systune),t3,'r','LineWidth',1.5)
plot(Ni,t2,'b','LineWidth',1.5)
plot(Ni(1:N_lmi),t4,'color',[0 0.5 0],'LineWidth',1.5)
legend('dlqr','systune','Alg. 1','LMI')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
ylim([0 60])
