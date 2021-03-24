function [K, f]=alg1_dof(P,S,nkx)
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
    N=size(S,1);
    nyi=ny/N;
    nui=nu/N;
    I=eye(N);
    A=P.A;
    B1=P.B(:,1:nd);
    B2=P.B(:,nd+1:end);
    C1=P.C(1:nz,:);
    C2=P.C(nz+1:end,:);
    D12=P.D(1:nz,nd+1:end);
    D21=P.D(nz+1:end,1:nd);
    NKX=N*nkx;
    
    At=[A zeros(nx,NKX); zeros(NKX,nx) zeros(NKX)];
    B1t=[B1; zeros(NKX,nd)];
    B2t=[zeros(nx,NKX) B2; eye(NKX) zeros(NKX,nu)];
    C1t=[C1 zeros(nz,NKX)];
    D12t=[zeros(nz,NKX) D12];
    C2t=[zeros(NKX,nx) eye(NKX); C2 zeros(ny,NKX)];
    D21t=[zeros(NKX,nd); D21];
    
    SAk=kron(S,ones(nkx));
    SBk=kron(I,ones(nkx,nyi));
    SCk=kron(S,ones(nui,nkx));
    SDk=kron(I,ones(nui,nyi));
    St=[SAk SBk; SCk SDk];
    
    pars.plantinfo.A=At;
    pars.plantinfo.B1=B1t;
    pars.plantinfo.B=B2t;
    pars.plantinfo.C1=C1t;
    pars.plantinfo.C=C2t;
    pars.plantinfo.D12=D12t;
    pars.plantinfo.D21=D21t;
    pars.plantinfo.S=St;
    
    % check if A is stable
    nvar=sum(sum(St));
    x0=0.1*randn(nvar,1);
    
    K0=zeros(nu+NKX,ny+NKX);
    count=1;
    for j=1:ny+NKX
        for i=1:nu+NKX
            if St(i,j)==1
                K0(i,j)=x0(count);
                count=count+1;
            end
        end
    end
    
    Acl=At+B2t*K0*C2t;
    specrad=max(abs(eig(Acl)));
    if specrad < 1
        stable=1;
    else
        stable=0;
    end 

    pars.nvar=nvar;
    options.nstart=1;
    options.x0=x0;
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
    Km=zeros(nu+NKX,ny+NKX);
    count=1;
    for j=1:ny+NKX
        for i=1:nu+NKX
            if St(i,j)==1
                Km(i,j)=x(count);
                count=count+1;
            end
        end
    end
    AK=Km(1:NKX,1:NKX);
    BK=Km(1:NKX,NKX+1:end);
    CK=Km(NKX+1:end,1:NKX);
    DK=Km(NKX+1:end,NKX+1:end);
    
    K=ss(AK,BK,CK,DK,-1);   
end



