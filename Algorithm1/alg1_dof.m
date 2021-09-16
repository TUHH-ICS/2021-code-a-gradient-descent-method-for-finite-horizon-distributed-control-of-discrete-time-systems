function [K] = alg1_dof(A,B,C,Q,QT,R,Vi,alpha0,Nmax,Tmax,S,nkx)
% Algorithm 1 with dynamic output feedback

% input:
%   - sparse system matrices A, B and C
%   - sparse cost function matrices Q, QT and R
%   - variance of measurement noise
%   - initial step length for gradient descent: alpha0
%   - maximum number of iterations: Nmax
%   - length of the finite time horizon: Tmax
%   - sparsity pattern of the state feedback matrix: S
%   - subsystem order of the controller

% output:
%   - sparse state feedback gain

% Assumptions: B,Q,R are blockdiagonal, Qxu=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx=size(A,1);
nu=size(B,2);
ny=size(C,1);
N=size(S,1);
nui=nu/N;
nyi=ny/N;

nk=N*nkx;
I=speye(N);
% generate sparsity patterns
SA=kron(S,ones(nkx,nkx));
SB=kron(I,ones(nkx,nyi));
SC=kron(S,ones(nui,nkx));
SD=kron(I,ones(nui,nyi));

% generate initial random controller
AK=0.1*SA.*randn(nk,nk);
BK=0.1*SB.*randn(nk,ny);
CK=0.1*SC.*randn(nu,nk);
DK=0.1*SD.*randn(nu,ny);

i=1;
while (i<=Nmax)
    alpha=alpha0;
    x=zeros(nx,Tmax);
    xc=zeros(nk,Tmax);
    u=zeros(nu,Tmax);
    y=zeros(ny,Tmax);
    ax=zeros(nx,Tmax-1);
    axk=zeros(nk,Tmax-1);
    v=normrnd(0,Vi,[ny, Tmax]);
    x(:,1)=normrnd(0,1,nx,1);
    y(:,1)=C*x(:,1)+v(:,1);
    xc(:,1)=zeros(nk,1);
    u(:,1)=CK*xc(:,1)+DK*y(:,1);
    
    % step 1: simulate x, xc, u, y
    V0=0;
    for j=1:Tmax-1
        x(:,j+1)=A*x(:,j)+B*u(:,j);
        y(:,j+1)=C*x(:,j+1)+v(:,j+1);
        xc(:,j+1)=AK*xc(:,j)+BK*y(:,j);
        u(:,j+1)=CK*xc(:,j+1)+DK*y(:,j+1);
        V0=V0+x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j);
    end
    V0=V0+x(:,Tmax)'*QT*x(:,Tmax);
    
    % step 2: simulate the adjoint state
    ax(:,Tmax-1)=-QT*x(:,Tmax);
    axk(:,Tmax-1)=zeros(nk,1);
    for j=Tmax-2:-1:1
        ax(:,j)=(A+B*DK*C)'*ax(:,j+1)-Q*x(:,j+1)-C'*DK'*R*u(:,j+1)+C'*BK'*axk(:,j+1);
        axk(:,j)=AK'*axk(:,j+1)-CK'*R*u(:,j+1)+CK'*B'*ax(:,j+1);
    end
    
    % step 3: calculate the gradient
    GA=SA;
    [ii,jj] = find(GA);
    tmp2=xc(:,1:Tmax-1)';
    for k=1:length(ii)
        GA(ii(k),jj(k))=-2*axk(ii(k),:)*tmp2(:,jj(k));
    end
    
    GB=SB;
    [ii,jj] = find(GB);
    tmp2=y(:,1:Tmax-1)';
    for k=1:length(ii)
        GB(ii(k),jj(k))=-2*axk(ii(k),:)*tmp2(:,jj(k));
    end
    
    GC=SC;
    [ii,jj] = find(GC);
    tmp1=2*(R*u(:,1:Tmax-1)-B'*ax);
    tmp2=xc(:,1:Tmax-1)';
    for k=1:length(ii)
        GC(ii(k),jj(k))=tmp1(ii(k),:)*tmp2(:,jj(k));
    end
    
    GD=SD;
    [ii,jj] = find(GD);
    tmp1=2*(R*u(:,1:Tmax-1)-B'*ax);
    tmp2=y(:,1:Tmax-1)';
    for k=1:length(ii)
        GD(ii(k),jj(k))=tmp1(ii(k),:)*tmp2(:,jj(k));
    end
    
    % backtracking line search
    stop=1;
    while stop
        AKt=AK-alpha*GA;
        BKt=BK-alpha*GB;
        CKt=CK-alpha*GC;
        DKt=DK-alpha*GD;
        % Simulate using updated controller
        u(:,1)=CKt*xc(:,1)+DKt*y(:,1);
        Vn=0;
        for j=1:Tmax-1
            x(:,j+1)=A*x(:,j)+B*u(:,j);
            y(:,j+1)=C*x(:,j+1)+v(:,j+1);
            xc(:,j+1)=AKt*xc(:,j)+BKt*y(:,j);
            u(:,j+1)=CKt*xc(:,j+1)+DKt*y(:,j+1);
            Vn=Vn+x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j);
        end
        Vn=Vn+x(:,Tmax)'*QT*x(:,Tmax);
        
        if Vn<V0
            i=i+1;
            AK=AK-alpha*GA;
            BK=BK-alpha*GB;
            CK=CK-alpha*GC;
            DK=DK-alpha*GD;
            stop=0;
        else
            alpha=alpha/10;
            if alpha<1e-20
                i=Nmax+1;
                stop=0;
            end
        end
    end
end

K=ss(AK,BK,CK,DK,-1);
end







