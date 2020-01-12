function A = computeProbAdjacency(positions, covariances, varargin)
% COMPUTEPROBADJACENCY Compute an adjacency matrix using the probabilistic
% mass of the cells' covariance ellipses.
%
% Input:
%   positions     = cell positions
%   covariances   = cell positional covariances
%   prob_ellipse  = probabilistic mass for the covariance ellipses
%                   default = 0.95 (95%)
%   sample_points = number of points for estimating covariance ellipses
%
% Output:
%   A = adjacency matrix

% What is the probablistic mass for the covariance ellipses?
prob_ellipse = 0.95;
if ~isempty(varargin)
    prob_ellipse = varargin{1};
end

% How many points should we use to estimate the covariance ellipses?
sample_points = 1000;
if length(varargin) > 1
    sample_points = varargin{2};
end

% Use a stable seed for the random number generator.
rng('default');
rng(1);

% Initialize the data.
A = zeros(size(positions,1)); % adjacency matrix A
points = cell(1,size(A,1)); % points estimating the covariance ellipses

% Generate probabilistic masses for the cells' covariance.
for i = 1:size(A,1)
    
    % Generate the covariance ellipse sampling points randomly.
    points{i} = mvnrnd(positions(i,:), covariances(:,:,i), ...
        sample_points);
    % Sort the ellipsoid distances from the cell's center.
    d = pdist2(points{i}, positions(i,:), 'mahalanobis', ...
        covariances(:,:,i));
    [~,idx] = sort(d, 'ascend');
    % Sample points for the probabilistic mass of the cell's covariance.
    points{i} = points{i}(idx(1:sample_points * prob_ellipse),:);
end

% Compute the distances between neurons' potential positions.
disp('Computing distances between neurons, please wait patiently ...');
for i = 1:size(A,1)
    
    % Compute the distance from cell's point cloud to its mean.
    di = pdist2(points{i}, positions(i,:), 'mahalanobis', ...
        covariances(:,:,i));
    % Compare cell point clouds, pairwise.
    for j = 1:size(A,2)
        dj = pdist2(points{i}, positions(j,:), 'mahalanobis', ...
            covariances(:,:,j));
        if any(dj < di, 1)
            A(i,j)=1;
        end
    end
end

% Restore the random number generator.
rng('shuffle');
end
