function [X]=optimal_color_solver(A,R,k,margin,sparsity,iter,varargin)
% Solver for fractional graph coloring with k colors and rank constraints
% Solves <A,squareform(pdist(R*X))> + gamma * ||X||_1 s.t. R*X*1 = 1, X>=0
% Input: A - nxn adjacency matrix (may or may not be weighted, 0 on diag )
%        R - Rank constraint for colors i.e. colors Z have to be Z = R*X
%        k - number of colors
%        sparsity - sparsity level of X
%        iter - number of multiplicative updates
% Optional input:
%        is_plot - are we plotting convergence to the objective?
%           default = true
%        is_normalize - compute R1=1 or R1<=1?
%           default = true (R1=1)
% Output: X - Color components such that R*X = Z --> Z are the colors.

% Plot convergence to the objective?
is_plot = true;
if ~isempty(varargin)
    is_plot = varargin{1};
end

% Compute R1=1 or R1<=1?
% If is_normalize=true, then colors sum to 1.
is_normalize = true;
if length(varargin) > 1
    is_normalize = varargin{2}; 
end

% Solve the objective.
lambda=128;
gamma=sparsity;
vec=@(x)(x(:));
X=rand(size(R,2),k);
Z=rand(size(A,1),k);
for t=1:iter
    for i=1:size(A,1)
        
        dxp(i,:)=zeros(1,k);
        dxn(i,:)=zeros(1,k);
        for j=1:size(A,2)
            if and(A(i,j)>0,norm(Z(i,:)-Z(j,:))<margin)
                dxp(i,:)=dxp(i,:)+2*A(i,j)*Z(i,:);
                dxn(i,:)=dxn(i,:)+2*A(i,j)*Z(j,:);
            end
        end
    end
    Z = (Z./(dxn+R*X+eps)).*(dxp + lambda*Z);
    Z=Z+eps;
    if is_normalize
        Z=Z./sum(Z,2);
    else
        Z=Z./(max(Z,[],2));
    end
    X = (X./(2*lambda*R'*R*X + gamma)).*(2*lambda*R'*Z);
    X=X./sum(X,1)*max(sum(X,1)); % Constraint that every column of X must be used to prevent column degeneration
    Z = R*X;
    
    % Plot convergence to the objective?
    if is_plot
        obj(t,1) = sum(vec(max(margin*A-squareform(pdist(R*X)).*A,0)));
        obj(t,2) = sum(vec(max(margin*A-squareform(pdist(R*X)).*A,0))>0);
        plot(obj)
        grid on
        legend({'Sum of margin violations','# of margin violators'});
        drawnow
    end
end
end
