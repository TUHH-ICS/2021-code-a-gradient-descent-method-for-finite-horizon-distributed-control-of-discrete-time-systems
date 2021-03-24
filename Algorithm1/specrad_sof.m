function [f, g] = specrad_sof(x,pars)
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
    C=plant.C;
    S=plant.S;
    nx=size(A,1);
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
        % fprintf('specabsc:****** matrix of right eigenvectors is singular, perturbing\n')
        V = V + perturb*randn(size(V));
        u = e/V;
        perturb = 10*perturb;
    end
    u = u'; % by convention the conjugate transpose is usually called left eigenvector
    % (right eigenvector of A' for the conjugate eigenvalue)
    % r = norm(A'*u-lam'*u); % these are very consistent in random tests
    % r_alt = norm(A'*u_alt-lam'*u_alt);  % these vary widely
    % G = real(conj(lam)/f * u*v');  % gradient in complex matrix space
    G = real( conj(lam)/abs(lam) * conj(u)*conj(v)');
    GK=B'*G*C';
    count=1;
    g=zeros(size(x));
    for j=1:ny
        for i=1:nu
            if S(i,j)==1
                g(count)=GK(i,j);
                count=count+1;
            end
        end
    end
    
end


