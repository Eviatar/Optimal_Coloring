function [Z]=optimal_color_solver(A,R,k,margin,sparsity,iter,varargin)
% Solver for fractional graph coloring with k colors and rank constraints
% Solves <A,squareform(pdist(R*Z))> + gamma * ||Z||_1 s.t. R*Z*1 = 1, Z>=0
% Input: A - nxn adjacency matrix (may or may not be weighted, 0 on diag )
%        R - Rank constraint for colors i.e. colors X have to be X = R*Z
%        k - number of colors
%        sparsity - sparsity level of Z
%        iter - number of multiplicative updates
% Output: Z - Color components such that R*Z = X --> X are the colors.

% Plot convergence to the objective?
is_plot = true;
if ~isempty(varargin)
    is_plot = varargin{1};
end

% Solve the objective.
lambda=128;
gamma=sparsity;
vec=@(x)(x(:));
X=rand(size(A,1),k);
Z=rand(size(R,2),k);
for t=1:iter
    for i=1:size(A,1)
        
        dxp(i,:)=zeros(1,k);
        dxn(i,:)=zeros(1,k);
        for j=1:size(A,2)
            if and(A(i,j)>0,norm(X(i,:)-X(j,:))<margin)
            dxp(i,:)=dxp(i,:)+2*A(i,j)*X(i,:);
            dxn(i,:)=dxn(i,:)+2*A(i,j)*X(j,:);
            end
        end
    end
    X = (X./(dxn+R*Z+eps)).*(dxp + lambda*X);
    X=X+eps;
    X=X./(sum(X,2));
    Z = (Z./(2*lambda*R'*R*Z + gamma)).*(2*lambda*R'*X);
    Z=Z./sum(Z,1)*max(sum(Z,1)); % Constraint that every column of Z must be used to prevent column degeneration
    X = R*Z;
    
    % Plot convergence to the objective?
    if is_plot
        obj(t,1) = sum(vec(max(margin*A-squareform(pdist(R*Z)).*A,0)));
        obj(t,2) = sum(vec(max(margin*A-squareform(pdist(R*Z)).*A,0))>0);
        plot(obj)
        grid on
        legend({'Sum of margin violations','# of margin violators'});
        drawnow
    end
end
end
