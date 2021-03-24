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
else
    N=150;
end

% dimensions of subsystems
nx=3;
nu=1;
% sparsity patterns
S=eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
Sa=kron(S,ones(nx));
Sb=kron(S,ones(nx,nu));

for i=1:NoS
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
    Ts=-1;
    gen_plant=ss(A,[eye(nx*N) B],[sqrt(Q); zeros(N,nx*N); eye(nx*N)],[zeros(nx*N) zeros(nx*N,N); zeros(N,nx*N) sqrt(R); zeros(nx*N,(nx+1)*N)],Ts);
    
    % dlqr:
    K1=dlqr(A,B,Q,R);
    cost1(i)=norm(lft(gen_plant,-K1));
    % Algorithm 2:
    Tmax=25; 
    Nmax=1000; 
    alpha0=0.01;
    
    % convert dense matrices to sparse matrices
    As=sparse(A);
    Bs=sparse(B);
    Qs=sparse(Q);
    QTs=sparse(QT);
    Rs=sparse(R);
    Ss=sparse(kron(S, ones(nu,nx)));
 
    K3=alg2_sf(As,Bs,Qs,QTs,Rs,Nmax,Tmax,Ss,alpha0);
    cost3(i)=norm(lft(gen_plant,-K3));
    if Case==1
        % Algorithm 1:
        [K2, f] = alg1_sf(A,B,Q,zeros(N*nx,N*nu),R,kron(S,ones(nu,nx)));
        cost2(i)=norm(lft(gen_plant,-K2));
        % upper bound using System Level Synthesis
        % first an FIR constraint of 20 and then an FIR constraint of 30 is choosen.
        % Since the maximum difference is less then 1e-3, we conclude
        % the SLS method converged
        Sys=LTISystem;
        Sys.Nx=N*nx;
        Sys.Nw=N*nx;
        Sys.Nu=N*nu;
        Sys.Nz=N*(nx+nu);
        Sys.A=As;
        Sys.B1=speye(N*nx);
        Sys.B2=Bs;
        Sys.C1=[sqrt(Qs); sparse(nu*N,nx*N)];
        Sys.D12=[sparse(N*nx,N*nu); sqrt(R)];
        Sys.D11=sparse(Sys.Nz,Sys.Nx);
        
        slsParams=SLSParams();
        slsParams.T_=30;
        slsParams.add_constraint(SLSConstraint.CommSpeed, 1);
        slsParams.add_objective(SLSObjective.H2, 1);
        clMaps=state_fdbk_sls(Sys,slsParams);
        Ksls = Ctrller.ctrller_from_cl_maps(clMaps);
        lqrCost = get_ctrller_lqr_cost(Sys, Ksls, 16);
        costSLS30(i)=sqrt(lqrCost);
        slsParams.T_ = 40;
        clMaps=state_fdbk_sls(Sys,slsParams);
        Ksls = Ctrller.ctrller_from_cl_maps(clMaps);
        lqrCost = get_ctrller_lqr_cost(Sys, Ksls, 16);
        costSLS40(i)=sqrt(lqrCost);
        
        % LMI method
        P.A=A;
        P.B1=eye(nx*N);
        P.B2=B;
        P.C1=[sqrt(Q); zeros(N,nx*N)];
        P.D12=[zeros(nx*N,N); sqrt(R)];

        K4 = LMI_sf(P,N,kron(S,ones(nu,nx)));
        cost4(i)=norm(lft(gen_plant,K4));
    end


       
end
%%
if Case==1
    disp('max difference between the cost of the SLS using a horizon of 30 and 40')
    max(abs(costSLS40-costSLS30))

    G2=cost1./cost2;
    G3=cost1./cost3;
    G4=cost1./cost4;
    G5=cost1./costSLS40;

    G2_=fliplr(sort([G2 1]));
    G3_=fliplr(sort([G3 1]));
    G4_=fliplr(sort([G4 1]));
    G5_=fliplr(sort([G5 1]));
else
    G=cost1./cost3;
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
    legend('Alg. 1','Alg. 2','LMI','upper bound')
    xlabel('Fraction of systems [%]')
    ylabel('Centralized/ distribtuted')
else
    figure(1)
    plot(100*x_,G_,'b','LineWidth',1.5)
    grid on
    box on
    legend('Alg. 2')
    xlabel('Fraction of systems [%]')
    ylabel('Centralized/ distribtuted')
    ylim([0 1])
end
