%% Clean up the workspace.
clear all;
clc;
close all;


%% Load the user parameters.
optimal_color_parameters;


%% Plot customization.

% Are we plotting the body as a whole or in parts?
is_plot_whole_body = true;

% Are we plotting the color margin violators?
is_plot_color_violators = false;


%% Helper functions.
vec=@(x)(x(:));


%%  Clean up the reporters.

% Remove all or none cell reporters.
all_cells = all(R == 1, 1);
R(:,all_cells) = [];
reporter_names = reporter_names(~all_cells);
no_cells = all(R == 0, 1);
R(:,no_cells) = [];
reporter_names = reporter_names(~no_cells);

% Remove redundant reporters.
[~,unique_R,~] = unique(double(R>0)','rows', 'stable');
R = R(:,unique_R);
reporter_names = reporter_names(unique_R);


%% KNN parameters -- NOT USED.
%knn = 2; % KNN graph neighbors
%Ak=zeros(size(cell_distances));
%for i=1:size(cell_distances,1)
%    [~,idx]=sort(cell_distances(i,:),'ascend');
%    Ak(i,idx(1:knn+1))=1;
%end
%Ak=Ak-diag(diag(Ak));
% A = Ak; % knn graph


%% Solve the max margin fractional graph coloring.
%rng('default'), rng(1); % use a stable seed for the random number generator
X=optimal_color_solver(A,R,num_colors,color_margin,sparsity,iterations);


%% Do we have the NeuroPAL info?
is_NeuroPAL = false;
if exist('NeuroPAL_colors', 'var')
    is_NeuroPAL = true;
end


%% Compute the real colors and reporters.

% Compute the optimal reporters.
%optimal_reporters = find(sum(X,2)>0.1); % use a hard threshold of 0.1
optimal_reporters = find(sum(X,2) > graythresh(X)); % use Otsu to choose the threshold

% Compute the optimal colors.
optimal_colors = R(:,optimal_reporters) * X(optimal_reporters,:);
optimal_colors = optimal_colors ./ max(optimal_colors);
optimal_colors(isnan(optimal_colors)) = 0; % 0 the NaNs

% Compute the NeuroPAL colors.
if is_NeuroPAL
    NeuroPAL_colors = NeuroPAL_colors ./ max(NeuroPAL_colors);
    NeuroPAL_colors(isnan(NeuroPAL_colors)) = 0; % 0 the NaNs
end


%% Group the neurons for visualization.

% Initialize the head ganglia.
head_ganglia = 1:9;
head_ganglia_left = [1,2,3,5,6,9];
head_ganglia_right = [1,2,4,5,7,9];
head_ganglia_ventral = [1,2,3,4,5,8,9];

% Initialize the tail ganglia.
tail_ganglia = 14:17;
tail_ganglia_left = [14,15,16];
tail_ganglia_right = [14,15,17];
tail_ganglia_ventral = [14,15,16,17];

% Extra neurons to show.
extra_head_neurons = {'VA2';'VB3';'AS2'};
extra_tail_neurons = {'AS10';'DA7';'VD11';'VA11'};

% Initialize the head neurons.
head_neurons = vertcat(ganglia(head_ganglia).neurons, extra_head_neurons);
[~,head_i,~] = intersect(neurons, head_neurons);
head_neurons_left = vertcat(ganglia(head_ganglia_left).neurons, ...
    extra_head_neurons);
[~,head_left_i,~] = intersect(neurons, head_neurons_left);
head_neurons_right = vertcat(ganglia(head_ganglia_right).neurons, ...
    extra_head_neurons);
[~,head_right_i,~] = intersect(neurons, head_neurons_right);
head_neurons_ventral = vertcat(ganglia(head_ganglia_ventral).neurons, ...
    extra_head_neurons);
[~,head_ventral_i,~] = intersect(neurons, head_neurons_ventral);

% Initialize the tail neurons.
tail_neurons = vertcat(ganglia(tail_ganglia).neurons, extra_tail_neurons);
[~,tail_i,~] = intersect(neurons, tail_neurons);
tail_neurons_left = vertcat(ganglia(tail_ganglia_left).neurons, ...
    extra_tail_neurons);
[~,tail_left_i,~] = intersect(neurons, tail_neurons_left);
tail_neurons_right = vertcat(ganglia(tail_ganglia_right).neurons, ...
    extra_tail_neurons);
[~,tail_right_i,~] = intersect(neurons, tail_neurons_right);
tail_neurons_ventral = vertcat(ganglia(tail_ganglia_ventral).neurons, ...
    extra_tail_neurons);
[~,tail_ventral_i,~] = intersect(neurons, tail_neurons_ventral);


%% Visualize the results.

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Arial');
%set(0,'DefaultAxesFontSize', 18);

% Change default text fonts.
set(0,'DefaultTextFontname', 'Arial');
%set(0,'DefaultTextFontSize', 18);


%% Plot the optimal color graphs.
figure('Name', 'Color Margin Violation Graphs', 'NumberTitle',  'off');

% Do we have the NeuroPAL info?
num_plots = 2;
if is_NeuroPAL
    num_plots = 3;
end

% Graph structure.
subplot(1,num_plots,1);
colormap('gray');
imagesc(A);
colorbar;
title('Graph structure');
set(gca, 'Xtick',(1:length(cell_names)), 'XTickLabel', cell_names, ...
    'YTick', (1:length(cell_names)), 'YTickLabel', cell_names, ...
    'XTickLabelRotation', 90);

% Optimal color margin violations.
subplot(1,num_plots,2);
optimal_color_graph = squareform(pdist(optimal_colors));
optimal_color_margin = sum(vec(max(color_margin*A-optimal_color_graph.*A,0)))./sum(vec(A));
imagesc(optimal_color_graph.*A);
colorbar;
set(gca, 'Xtick', (1:length(cell_names)), 'XTickLabel', cell_names, ...
    'YTick', (1:length(cell_names)), 'YTickLabel', cell_names, ...
    'XTickLabelRotation', 90);
title(['Approximately optimal color margin violation: ' ...
    num2str(optimal_color_margin * 100, '%2.1f') '%']);

% NeuroPAL color margin violations.
if is_NeuroPAL
    subplot(1,num_plots,3);
    neuropal_color_graph = squareform(pdist(NeuroPAL_colors));
    neuropal_color_margin = sum(vec(max(color_margin*A-neuropal_color_graph.*A,0)))./sum(vec(A));
    imagesc(neuropal_color_graph.*A);
    colorbar;
    set(gca, 'Xtick', (1:length(cell_names)), 'XTickLabel', cell_names, ...
        'YTick', (1:length(cell_names)), 'YTickLabel', cell_names, ...
        'XTickLabelRotation', 90);
    title(['NeuroPAL color margin violation: ' ...
        num2str(neuropal_color_margin * 100, '%2.1f') '%']);
end


%% Plot the optimal vs. NeuroPAL colors.

% Plot the optimal colors.
fig_title = 'Near Optimal Colors';
plot_title = ['Near Optimal Colors - color margin violation = ' ...
    num2str(optimal_color_margin * 100, '%2.1f') '%'];

% Plot the whole body.
if is_plot_whole_body
plotCells(neurons, neuron_positions, optimal_colors, [], [], [], ...
    fig_title, plot_title);
else
    
% Plot the head.
plotCells(neurons(head_left_i), neuron_positions(head_left_i,:), ...
    optimal_colors(head_left_i,:), [], [], [0,-180], ...
    ['Head (Left) - ' fig_title], ['Head (Left) - ' plot_title]);
plotCells(neurons(head_right_i), neuron_positions(head_right_i,:), ...
    optimal_colors(head_right_i,:), [], [], [0,-180], ...
    ['Head (Right) - ' fig_title], ['Head (Right) - ' plot_title]);
plotCells(neurons(head_ventral_i), neuron_positions(head_ventral_i,:), ...
    optimal_colors(head_ventral_i,:), [], [], [], ...
    ['Head (Ventral) - ' fig_title], ['Head (Ventral) - ' plot_title]);

% Plot the tail.
plotCells(neurons(tail_left_i), neuron_positions(tail_left_i,:), ...
    optimal_colors(tail_left_i,:), [], [], [0,-180], ...
    ['Tail (Left) - ' fig_title], ['Tail (Left) - ' plot_title]);
plotCells(neurons(tail_right_i), neuron_positions(tail_right_i,:), ...
    optimal_colors(tail_right_i,:), [], [], [0,-180], ...
    ['Tail (Right) - ' fig_title], ['Tail (Right) - ' plot_title]);
plotCells(neurons(tail_ventral_i), neuron_positions(tail_ventral_i,:), ...
    optimal_colors(tail_ventral_i,:), [], [], [], ...
    ['Tail (Ventral) - ' fig_title], ['Tail (Ventral) - ' plot_title]);
end

% Plot the NeuroPAL colors.
if is_NeuroPAL
fig_title = 'NeuroPAL Colors';
plot_title = ['NeuroPAL Colors - color margin violation = ' ...
    num2str(neuropal_color_margin * 100, '%2.1f') '%'];

% Plot the whole body.
if is_plot_whole_body
plotCells(neurons, neuron_positions, NeuroPAL_colors, [], [], [], ...
    fig_title, plot_title);
else
    
% Plot the head.
plotCells(neurons(head_left_i), neuron_positions(head_left_i,:), ...
    NeuroPAL_colors(head_left_i,:), [], [], [0,-180], ...
    ['Head (Left) - ' fig_title], ['Head (Left) - ' plot_title]);
plotCells(neurons(head_right_i), neuron_positions(head_right_i,:), ...
    NeuroPAL_colors(head_right_i,:), [], [], [0,-180], ...
    ['Head (Right) - ' fig_title], ['Head (Right) - ' plot_title]);
plotCells(neurons(head_ventral_i), neuron_positions(head_ventral_i,:), ...
    NeuroPAL_colors(head_ventral_i,:), [], [], [], ...
    ['Head (Ventral) - ' fig_title], ['Head (Ventral) - ' plot_title]);

% Plot the tail.
plotCells(neurons(tail_left_i), neuron_positions(tail_left_i,:), ...
    NeuroPAL_colors(tail_left_i,:), [], [], [0,-180], ...
    ['Tail (Left) - ' fig_title], ['Tail (Left) - ' plot_title]);
plotCells(neurons(tail_right_i), neuron_positions(tail_right_i,:), ...
    NeuroPAL_colors(tail_right_i,:), [], [], [0,-180], ...
    ['Tail (Right) - ' fig_title], ['Tail (Right) - ' plot_title]);
plotCells(neurons(tail_ventral_i), neuron_positions(tail_ventral_i,:), ...
    NeuroPAL_colors(tail_ventral_i,:), [], [], [], ...
    ['Tail (Ventral) - ' fig_title], ['Tail (Ventral) - ' plot_title]);
end
end


%% Plot the optimal vs. NeuroPAL margin violators.
if is_plot_color_violators
% Plot the optimal color margin violators.
fig_title = 'Near Optimal Color Margin Violators';
plot_title = 'Near Optimal Colors';
plotCells(neurons, neuron_positions, optimal_colors, A, color_margin, ...
    [], fig_title, plot_title);

if is_NeuroPAL
fig_title = 'NeuroPAL Color Margin Violators';
plot_title = 'NeuroPAL Optimal Colors';
plotCells(neurons, neuron_positions, NeuroPAL_colors, A, color_margin, ...
    [], fig_title, plot_title);
end
end

%% Plot the reporters and their color assignments.
figure('Name', 'Optimal Reporters and Coloring', 'NumberTitle',  'off');
RGB_labels = {'Red','Green','Blue'};
colormap('gray');
imagesc(X(optimal_reporters,:));
c = colorbar;
c.Label.String = 'Reporter Concentration (%)';
c_lim = c.Limits;
c.Ticks = linspace(c_lim(1), c_lim(2), 11);
c.TickLabels = arrayfun(@(x) num2str(x), 0:10:100, 'UniformOutput', false);
set(gca,'Ytick', (1:length(optimal_reporters)), ...
    'Yticklabel', reporter_names(optimal_reporters), ...
    'Xtick', (1:num_colors));
if num_colors <= 3
    set(gca, 'Xticklabel', RGB_labels(1:num_colors));
end
xlabel('Reporter Colors');
title(['Reporters and Colors (Reporters = ' num2str(length(optimal_reporters)) ')']);
