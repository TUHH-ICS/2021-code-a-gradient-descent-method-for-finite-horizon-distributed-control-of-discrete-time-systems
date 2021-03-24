function [K, f] = alg1_sf(A,B,Qx,Qxu,Qu,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nx=size(A,1);
    nu=size(B,2);
    % check if A is stable
    specrad=max(abs(eig(A)));
    if specrad < 1
        stable=1;
    else
        stable=0;
    end 
    pars.plantinfo.A=A;
    pars.plantinfo.B=B;
    pars.plantinfo.Qx=Qx;
    pars.plantinfo.Qxu=Qxu;
    pars.plantinfo.Qu=Qu;
    pars.plantinfo.S=S;
    nvar=sum(sum(S));
    pars.nvar=nvar;
    options.nstart=1;
    options.x0=zeros(nvar,1);
    options.fvalquit=0.95;
    pars.fgname='specrad';
    % step 1: minimize spectral radius
    if ~stable
        
        [x, f]=bfgs(pars,options);
        
        if f>=1
            disp('unable to find stabilizing controller')
            K=[];
            f=inf;
            return
        end
    end
    % step 2: minimize LQR cost
    pars.fgname='DLQR_cost';
    options.x0=x;
    options.fvalquit=-inf;
    [x, f]=bfgs(pars,options);
    % construct controller
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
end

