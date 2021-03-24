function [f, g] = specrad(x,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
%
%  modified version based on
%
%  HIFOO, A Matlab package for Fixed Order H-infinity control
%  Copyright (C) 2008  Marc Millstone, Michael Overton
%
%  Discrete-time Modifications: 2009, Andrey Popov
%
%  Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plant=pars.plantinfo;
    A=plant.A;
    B=plant.B;
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
    [V, Lambda] = eig(Acl);
    lambda = diag(Lambda);
    [f, k] = max(abs(lambda));
    lam = lambda(k); % corresponding eigenvalue
    v = V(:,k);  % and right eigenvector
    n = length(A);
    I = eye(n);
    e = I(k,:);
    u = e/V;   % relevant row of inverse of V ("left row eigenvector")
    % might possibly have to perturb V to make it nonsingular
    perturb = 1e-16*max(max(abs(V)));
    while isnan(u)
        V = V + perturb*randn(size(V));
        u = e/V;
        perturb = 10*perturb;
    end
    u = u'; 
    G = real( conj(lam)/abs(lam) * conj(u)*conj(v)');
    GK=-B'*G;
    count=1;
    g=zeros(size(x));
    for j=1:nx
        for i=1:nu
            if S(i,j)==1
                g(count)=GK(i,j);
                count=count+1;
            end
        end
    end
    
end

