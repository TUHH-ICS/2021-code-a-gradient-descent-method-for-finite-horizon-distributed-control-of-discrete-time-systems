function K = alg1_sf(A,B,Q,QT,R,Nmax,Tmax,S,alpha0,K0)
% Algorithm 1 with state feedback

% input:
%   - sparse system matrices A and B
%   - sparse cost function matrices Q, QT and R
%   - initial step length for gradient descent: alpha0
%   - maximum number of iterations: Nmax
%   - length of the finite time horizon: Tmax
%   - sparsity pattern of the state feedback matrix: S
%   - (optional) initial controller (default: K=0)

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

if exist('K0','var')
    K=K0;
else
    K=sparse(nu,nx);
end

i=1;
while (i<=Nmax)
    alpha=alpha0;
    x=zeros(nx,Tmax);
    u=zeros(nu,Tmax);
    y=zeros(nx,Tmax-1);
    x(:,1)=normrnd(0,1,nx,1);
    x(:,1)=x(:,1);
    u(:,1)=K*x(:,1);
    Vo=0;
    % step 1: simulate x,u
    for j=1:Tmax-1
        x(:,j+1)=A*x(:,j)+B*u(:,j);
        u(:,j+1)=K*x(:,j+1);
        Vo=Vo+x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j);
    end
    Vo=Vo+x(:,Tmax)'*QT*x(:,Tmax);
    
    % step 2: simulate the adjoint state
    y(:,Tmax-1)=-QT*x(:,Tmax);
    for j=Tmax-2:-1:1
        y(:,j)=(A+B*K)'*y(:,j+1)-Q*x(:,j+1)-K'*R*u(:,j+1);
    end
    
    % step 3: calculate the gradient
    tmp=2*(R*u(:,1:Tmax-1)-B'*y);
    tmp2=x(:,1:Tmax-1)';
    G=S;
    [ii,jj] = find(G);

    for k=1:length(ii)
        G(ii(k),jj(k))=tmp(ii(k),:)*tmp2(:,jj(k));
    end
    
    % backtracking line search
    stop=1;
    while stop
        Ktmp=K-alpha*G;
        % Simulate using updated K
        u(:,1)=Ktmp*x(:,1);
        Vn=0;

        for j=1:Tmax-1
            x(:,j+1)=A*x(:,j)+B*u(:,j);
            u(:,j+1)=Ktmp*x(:,j+1);
            Vn=Vn+x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j);
        end
        Vn=Vn+x(:,Tmax)'*QT*x(:,Tmax);
        if Vn<Vo
            i=i+1;
            K=K-alpha*G;
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


end





