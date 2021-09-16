function K = LMI_sf(P,N,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Uses an own implementaion of the work of 
% A. Vamsi and N. Elia, "Design of Distributed Controllers realizable over arbitrary networks", 
% IEEE Conference on Decision and Control, 2010.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    A=P.A;
    B1=P.B1;
    B2=P.B2;
    C1=P.C1;
    D12=P.D12;
    
    nx=size(A,1);
    nu=size(B2,2);
    nz=size(C1,1);
    nd=size(B1,2);
    
    nxi=nx/N;
    
    P=sdpvar(nx);
    W=sdpvar(nz);
    X=[];
    for i=1:N
        X=blkdiag(X,sdpvar(nxi,nxi,'full'));
    end
    L=sdpvar(nu,nx,'full');
    for i=1:nx
        for j=1:nu
            if S(j,i)==0
                L(j,i)=0;
            end
        end
    end
    
    L1=[     P           A*X+B2*L     B1; ...
        (A*X+B2*L)'      X+X'-P     zeros(nx,nd); ...
            B1'      zeros(nd,nx)    eye(nd)];
        
    L2=[W  C1*X+D12*L; (C1*X+D12*L)'  X+X'-P];
    
    LMI=[L1>0, L2>0];
    obj=trace(W);
    
    OPTIONS = sdpsettings;
    OPTIONS.solver='sedumi';
    diagnostics = optimize(LMI,obj,OPTIONS);
    
    K=value(L)/value(X);         
end

