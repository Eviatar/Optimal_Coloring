%% PLEASE READ ALL OF ME!!!

%% FILE DESCRIPTION.
% This file contains the parameters for the script "find_optimal_colors.m".
% The accompanying NeuroPAL paper has full details for the algorithm.
% The algorithm chooses a near-optimal set of reporters and assigns them
% colors in order to provide an unique color barcode that discriminates
% individual cell types.
%
% Provided a sufficiently large dataset, the algorithm produces a different
% near-optimal solution each time it is run. These solutions can be tested
% in vivo to choose the best one and/or to help refine the algorithmic
% parameters using the real-world results.

%% EXAMPLE PARAMETERS.
% We have included an example application below: a near-optimal NeuroPAL.
% To run the example:
% 1. Place all accompanying files in a single directory.
% 2. Open Matlab.
% 3. Navigate to the directory with the files.
% 4. Execute the script "find_optimal_colors.m"

%% USER PARAMETERS.
% Please adapt the following parameters for your needs. Each parameter has
% an accompanying decription to guide you in choosing appropriate values.

%% Load the Brain Atlas and NeuroPAL-tested reporter data.
% This dataset includes:
% 1. BrainAtlas = a table of the Brain Atlas reporters
% 2. NeuroPAL = a table of the NeuroPAL-tested reporters
% 3. neurons = the names of the head and tail neurons to color
% 4. neuron_distances = the distances between the neurons
% 5. neuron_positions = the physical positions for the neurons
% 6. NeuroPAL_colors = the NeuroPAL colors for the neurons
load data.mat;

%% PLEASE SET THE FOLLOWING USER PARAMETERS FOR YOUR APPLICATION.

%% The reporter matrix (R) of reporters x cellular expression.
%% NOTE: YOU MUST INCLUDE THE REPORTER AND CELL NAMES IN CELL ARRAYS!
% reporter_names = the reporter names
% cell_names = the cell (e.g., neuron) names
%
% A matrix with the rows = cells and columns = reporter expression.
% The matrix values can be binary: expressed = 0 and not expressed = 1.
% Or, the matrix can reflect the reporter weight per cell, normalized to
% [0,1] (e.g., one can normalize the transcript counts from scRNA-seq data).
%
% BrainAtlas uses binary weights for neural expression.
% NeuroPAL uses 5 quantized weights (0, 0.25, 0.5, 0.75, 1) for expression.
% Both tables are organized with the first 3 columns describing the neuron
% classes, names, and types, followed by the reporters and their
% expression weight per neuron.

% Choose the reporters.
reporter_table = NeuroPAL; % use the NeuroPAL-tested reporters
reporter_table = BrainAtlas; % use the Brain Atlas reporters
[reporter_names, cell_names, R] = getReporterInfo(reporter_table, neurons);

%% The adjacency matrix (A) for cells.
% A matrix of pairwise neurons that should be distinguishable.
% If 2 neurons should be distinguishable, set their matrix value in A = 1.
% Otherwise, set the matrix value in A = 0.

% NEIGHBORS = NEURONS WITH CENTERS WITHIN RADIAL DISTANCE
% Here we compute A using the distances between neurons. We indicate that
% neurons with 8 microns distance (~2 nuclei), should be distinguishable.
%radius = 8;
%A = double(neuron_distances < radius);
%A = A - diag(diag(A)); % a neuron should NOT be distinguishable from itself

% NEIGHBORS = NEURONS WHOSE PROBABILSITIC COVARIANCE ELLIPSE OVERLAPS
% Here we compute A using the probabilistic mass of the neurons' covariance
% ellipses. We indicate that neurons whose positional neighborhood has 95%
% probability of including each other, should be distinguishable.
if exist('neuron_prob_distances','var')
    A = neuron_prob_distances;

% Compute A using the probabilistic mass of the neurons' covariance ellipses.
else
    A = computeProbAdjacency(neuron_positions, neuron_covariances);
end

%% The number of colors (c) available.
% NeuroPAL uses 3 colors for its RGB pallette:
% 1. Red = mNeptune2.5
% 2. Green = CyOFP1
% 3. Blue = mTagBFP2
% TagRFP-T was used for panneuronal coloring but, you can repurpose it as a
% 4th color. Similarly, GFP was left empty to identify reporter signal but,
% you can reclaim it here as a 5th color.
% Note: we cannot visualize more than 3 colors on the screen.
num_colors = 3;

%% The color margin (intensity) to discriminate.
% If 2 neighboring neurons are 1 color margin's RGB distance away, then
% they are considered to be distinguishably colored.
% NeuroPAL reporters are stable enough to discriminate 3 levels of
% intensity (bright, medium, weak). Therefore, we use 1/3 as the margin. 
color_margin = 1/3; % bright, medium, weak

%% The sparsity for reporter coloring (lambda).
% High sparsity chooses fewer reporters for coloring. Conversely, low
% sparsity chooses more reporters for coloring.
% To match the NeuroPAL reporters set (40 RGB reporters + 1 panneuronal),
% we adjust the sparsity till the algorithm chooses ~40 reporters, lamdba =
% ~48. To limit the number of reporters (e.g., to choose ~3 reporters),
% increase the sparsity until the algorithm adresses your needs.
% For Brain Atlas: 50=~50, 100=~30 250=~15-25 500=~8-12 1000=~3-6 reporters
% For NeuroPAL: 0=~40, 50=~30
sparsity = 100; 

%% The number of algorithmic iterations to reach convergence.
% Running the algorithm plots our objective: minimizing the sum of color
% margin violations and the number of cells violating the color margin.
% In other words, minimizing the number of neighboring cells that are
% indistinguishably colored, while respecting the sparsity of reporters.
% The graph will rapidly decrease and eventually converge within some
% bounds. We suggest you start with 1000 iterations and adjust this to
% reflect a number of iterations beyond where you see convergence.
% For the Brain Atlas, convergence requires ~250 iterations.
% For the NeuroPAL-tested reporters, convergence requires ~100 iterations.
iterations = 300; % number of iterations in optimization, choose higher if not converged
