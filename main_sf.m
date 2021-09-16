%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
% Number of systems (can be decreased in order to reduce the computation time) 
NoS=50;
% Case=1: small systems consisting of 10 subsystems
% Case=2: large scale example consisting of 150 subsystems
Case=1;
if Case==1
    N=10;
    alp=0.3;
else
    N=150;
    alp=5/N;
end

% dimensions of subsystems
nx=3;
nu=1;

for i=1:NoS
    flag=1;
    % generate random sparsity patterns
    while(flag)
        St = rand(N) < alp;
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
    % rescaling of the A matrix
    A=1.5*A/max(abs(eig(A))); 
    % cost function
    Q=eye(nx*N);
    QT=0*Q;
    R=0.1*eye(nu*N);
    
    % generalized plant
    Ts=1;
    gen_plant=ss(A,[eye(nx*N) B],[sqrt(Q); zeros(N,nx*N); eye(nx*N)],[zeros(nx*N) zeros(nx*N,N); zeros(N,nx*N) sqrt(R); zeros(nx*N,(nx+1)*N)],Ts);

    % dlqr:
    K1=dlqr(A,B,Q,R);
    cost1(i)=norm(lft(gen_plant,-K1));
    % Algorithm 1:
    Tmax=25; 
    Nmax=1500; 
    alpha0=0.01;
    
    % convert dense matrices to sparse matrices
    As=sparse(A);
    Bs=sparse(B);
    Qs=sparse(Q);
    QTs=sparse(QT);
    Rs=sparse(R);
    Ss=sparse(kron(S, ones(nu,nx)));
    S1=sparse(S);
    Nx=nx*ones(1,N);
    Nu=nu*ones(1,N);
     
    K2=alg1_sf(As,Bs,Qs,QTs,Rs,Nmax,Tmax,Ss,alpha0);
    cost2(i)=norm(lft(gen_plant,K2))
    if Case==1
        % systune:
        K3=realp('K3',kron(S, ones(nu,nx)));
        K3.Free=kron(S, ones(nu,nx));
        CL0=lft(gen_plant,K3);
        CL0.InputName='in';
        CL0.OutputName='out';
        goal=TuningGoal.Variance('in','out',0.01);
        [CL,fSoft] = systune(CL0,goal);
        K3=CL.A.Blocks.K3.value;
        cost3(i)=norm(lft(gen_plant,K3));
        % SLS:
        [K4, obj]=SLS_sf(A,B,Q,R,double(S),25);
    	costSLS25(i)=norm(lft(gen_plant,K4),2);
        
        [K4, obj]=SLS_sf(A,B,Q,R,double(S),35);
    	costSLS35(i)=norm(lft(gen_plant,K4),2);
        
        %% LMI method:
        P.A=A;
        P.B1=eye(nx*N);
        P.B2=B;
        P.C1=[sqrt(Q); zeros(N,nx*N)];
        P.D12=[zeros(nx*N,N); sqrt(R)];

        K4 = LMI_sf(P,N,kron(S,ones(nu,nx)));
        cost5(i)=norm(lft(gen_plant,K4));
        
    end
end


%%
if Case==1
    disp('max difference between the cost of the SLS using a horizon of 35 and 25')
    max(abs((costSLS35-costSLS25)/costSLS35))

    G2=cost1./cost3;
    G3=cost1./cost2;
    G4=cost1./cost5;
    G5=cost1./costSLS35;

    G2_=fliplr(sort([G2 1]));
    G3_=fliplr(sort([G3 1]));
    G4_=fliplr(sort([G4 1]));
    G5_=fliplr(sort([G5 1]));
else
    G=cost1./cost2;
    G_=fliplr(sort([G 1]));
end
x_=0:(1/NoS):1;

if Case==1
    figure(1)
    plot(100*x_,G2_,'r','LineWidth',1.5)
    hold on
    grid on
    box on
    plot(100*x_,G3_,'b','LineWidth',1.5)
    plot(100*x_,G4_,'color',[0,0.5,0],'LineWidth',1.5)
    plot(100*x_,G5_,'--','color',[0.3,0.3,0.3],'LineWidth',1.5)
    legend('systune','Alg. 1','LMI','upper bound')
    xlabel('Fraction of systems [%]')
    ylabel('Centralized/ distribtuted')
    ylim([0 1])
else
    figure(1)
    plot(100*x_,G_,'b','LineWidth',1.5)
    grid on
    box on
    legend('Alg. 1')
    xlabel('Fraction of systems [%]')
    ylabel('Centralized/ distribtuted')
    ylim([0 1])
end
