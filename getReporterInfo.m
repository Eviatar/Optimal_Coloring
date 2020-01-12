function [reporters, neurons, expression] = ...
    getReporterInfo(reporter_table, varargin)
%GETREPORTERINFO Get the reporter information from a table.
%
% [REPORTERS, NEURONS, EXPRESSION] = GETREPORTERINFO(REPORTER_TABLE)
% [REPORTERS, NEURONS, EXPRESSION] = GETREPORTERINFO(REPORTER_TABLE,
%                                       NEURONS_TO_USE)
% [REPORTERS, NEURONS, EXPRESSION] = GETREPORTERINFO(REPORTER_TABLE,
%                                       NEURONS_TO_USE, IS_CLEAN_TABLE)
%
% Input:
%   reporter_table - the reporter table (BrainAtlas or NeuroPAL)
%   neurons_to_use - which neurons are we using?
%                    default: [] = all neurons
%   is_clean_table - are we cleaning up the table to remove redundant
%                    entries and entries with all or none expression?
%                    default: true
%
% Output:
%   reporters  - the reporter names for the expression columns
%   neurons    - the neuron names for the expression rows
%   expression - the expression matrix (neurons x reporters)

% Which neurons are we using?
neurons_to_use = [];
if ~isempty(varargin)
    neurons_to_use = varargin{1};
end

% Are we cleaning up the table?
is_clean_table = true;
if length(varargin) > 1
    is_clean_table = varargin{2};
end

% Set the reporter and neuron names.
reporters = reporter_table.Properties.VariableNames(4:end);

% Get the neurons we're using.
neurons = reporter_table.Neuron;
if isempty(neurons_to_use)
    
    % Extract the expression matrix.
    expression = double(reporter_table{:,4:end});
else
    
    % Extract the neurons we're using.
    use_neurons = categorical(neurons_to_use);
    [~, ~, n_i]=intersect(use_neurons, neurons, 'stable');
    neurons = reporter_table.Neuron(n_i);
    if ~all(use_neurons == neurons)
        error('Reporter weight matrix does NOT match the neuron list!');
    end
    
    % Extract the expression matrix.
    expression = double(reporter_table{n_i,4:end});
end

% Are we cleaning up the table?
if is_clean_table
    
    % Remove all or none neuron reporters.
    all_neurons = all(expression == 1, 1);
    expression(:,all_neurons) = [];
    reporters = reporters(~all_neurons);
    no_neurons = all(expression == 0, 1);
    expression(:,no_neurons) = [];
    reporters = reporters(~no_neurons);
    
    % Remove redundant reporters.
    [~,unique_reporters,~] = unique(double(expression>0)', 'rows', 'stable');
    expression = expression(:,unique_reporters);
    reporters = reporters(unique_reporters);
end
end
