function [K, f] = alg1_sof(P,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nx=size(P.A,1);
    nu=size(S,1);
    ny=size(S,2);
    nd=size(P.B,2)-nu;
    nz=size(P.C,1)-ny;
    
    pars.plantinfo.A=P.A;
    pars.plantinfo.B1=P.B(:,1:nd);
    pars.plantinfo.B=P.B(:,nd+1:end);
    pars.plantinfo.C1=P.C(1:nz,:);
    pars.plantinfo.C=P.C(nz+1:end,:);
    pars.plantinfo.D12=P.D(1:nz,nd+1:end);
    pars.plantinfo.D21=P.D(nz+1:end,1:nd);
    pars.plantinfo.S=S;
    
    % check if A is stable
    A=P.A;
    specrad=max(abs(eig(A)));
    if specrad < 1
        stable=1;
    else
        stable=0;
    end 
    nvar=sum(sum(S));
    pars.nvar=nvar;
    options.nstart=1;
    options.x0=zeros(nvar,1);
    options.fvalquit=0.95;
    pars.fgname='specrad_sof';
    % step 1: minimize spectral radius
    if ~stable
        [x, f]=bfgs(pars,options);
        options.x0=x;
        if f>=1
            disp('unable to find stabilizing controller')
            K=[];
            f=inf;
            return
        end
    end
    % step 2: minimize LQR cost
    pars.fgname='DH2_sof';
    options.fvalquit=-inf;
    [x, f]=bfgs(pars,options);
    
    % construct controller
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
    
end

