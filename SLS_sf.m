function [K, cost] = SLS_sf(A,B,Q,Ru,S,T)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
    % Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
    % Licensed under GPLv3. See License.txt in the project root for license information.
    % Uses an own implementaion of the work of 
    % J. Anderson, J. C. Doyle, S. H. Low and N. Matni, "System Level Synthesis", Annual Reviews in Control, vol. 47, pp. 364-393, 2019. 
    % Author: Simon Heinke
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % find support
    Ss=S;
    SnumM=0;
    for k=1:T
        St=(kron(Ss,ones(1,3))~=0);
        Msupp{k}=find(St);
        Ss=double(Ss*S>0);
        numM{k}=sum(sum(double(St)));
        SnumM=SnumM+numM{k};
    end
    
    N=size(B,2);
    Rsupp=find(ones(3*N,3*N));
    
    numR=9*N^2;
    count=SnumM+T*numR;

    spot=0;
    cvx_begin
        variable X(count)
        expression Rs(3*N,3*N,T)
        expression Ms(N,3*N,T)

        for t=1:T
            R{t}=Rs(:,:,t);
            R{t}(Rsupp)=X(spot+1:spot+numR);
            spot=spot+numR;

            M{t}=Ms(:,:,t);
            M{t}(Msupp{t})=X(spot+1:spot+numM{t});
            spot=spot+numM{t};
        end

        obj=0;
        for t=1:T
            vect=vec([sqrt(Q)*R{t}; sqrt(Ru)*M{t}]);
            obj=obj+vect'*vect;
        end

        minimize(obj)

        subject to
        R{1}==eye(3*N);
        for t=1:T-1
            R{t+1}==A*R{t}+B*M{t};
        end
        R{T}==zeros(3*N);
    cvx_end

    cost=sqrt(obj)
    
    %% controller construction
    Z=kron(diag(ones(T-2,1),-1),eye(N*3));
    I=zeros(N*3*(T-1),N*3);
    I(1:N*3,:)=eye(N*3);
    Rhat=[];
    Mhat=[];
    for t=1:T-1
        Rhat=[Rhat R{t+1}];
        Mhat=[Mhat M{t+1}];
    end
    M1=M{1};
    AK=Z-I*Rhat;
    BK=-I;
    CK=M1*Rhat-Mhat;
    DK=M1;

    K=ss(AK,BK,CK,DK,-1);
end

