function [f, g] = DLQR_cost(x,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plant=pars.plantinfo;
    A=plant.A;
    B=plant.B;
    Qx=plant.Qx;
    Qu=plant.Qu;
    Qxu=plant.Qxu;
    S=plant.S;
    nx=size(A,1);
    nu=size(B,2);
    
    K=zeros(nu,nx);
    count=1;
    for j=1:nx
        for i=1:nu
            if S(i,j)==1
                K(i,j)=x(count);
                count=count+1;
            end
        end
    end
    Acl=A-B*K;
    if max(abs(eig(Acl))) >=1
        f=inf;
        g=nan*ones(size(x));
        return
    end
    X0=dlyap(Acl,eye(nx));
    M=Qx-Qxu*K-K'*Qxu'+K'*Qu*K;
    P=dlyap(Acl',M);
    
    G=2*(Qu*K-Qxu'-B'*P*Acl)*X0;
    g=zeros(size(x));
    count=1;
    for j=1:nx
        for i=1:nu
            if S(i,j)==1
                g(count)=G(i,j);
                count=count+1;
            end
        end
    end
    
    f=trace((Qx-2*Qxu*K+K'*Qu*K)*X0);
end

