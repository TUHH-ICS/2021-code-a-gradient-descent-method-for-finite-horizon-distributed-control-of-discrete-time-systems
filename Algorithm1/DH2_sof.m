function [f, g] = DH2_sof(x,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plant=pars.plantinfo;
    A=plant.A;
    B1=plant.B1;
    B=plant.B;
    C1=plant.C1;
    C=plant.C;
    D12=plant.D12;
    D21=plant.D21;
    S=plant.S;
    
    nu=size(B,2);
    ny=size(C,1);
    
    K=zeros(nu,ny);
    count=1;
    for j=1:ny
        for i=1:nu
            if S(i,j)==1
                K(i,j)=x(count);
                count=count+1;
            end
        end
    end
    Acl=A+B*K*C;
    if max(abs(eig(Acl))) >=1
        f=inf;
        g=nan*ones(size(x));
        return
    end
    Bcl=B1+B*K*D21;
    Ccl=C1+D12*K*C;
    Dcl=D12*K*D21;
    
    Lc=dlyap(Acl,Bcl*Bcl');
    Lo=dlyap(Acl',Ccl'*Ccl);
    
    G=2*(B'*Lo*Acl*Lc*C'+B'*Lo*Bcl*D21'+D12'*Ccl*Lc*C'+D12'*Dcl*D21');
    
    g=zeros(size(x));
    count=1;
    for j=1:ny
        for i=1:nu
            if S(i,j)==1
                g(count)=G(i,j);
                count=count+1;
            end
        end
    end
    
    f=trace(Ccl*Lc*Ccl'+Dcl*Dcl');
end



