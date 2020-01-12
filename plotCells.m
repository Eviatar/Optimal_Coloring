function plotCells(names, positions, colors, varargin)
%PLOTCELLS Plot the cells in their colors.
%
% PLOTCELLS(NAMES, POSITIONS, COLORS)
% PLOTCELLS(NAMES, POSITIONS, COLORS, A, COLOR_MARGIN)
% PLOTCELLS(NAMES, POSITIONS, COLORS, A, COLOR_MARGIN, VIEW_ANGLE)
% PLOTCELLS(NAMES, POSITIONS, COLORS, A, COLOR_MARGIN, VIEW_ANGLE,
%           FIG_TITLE)
% PLOTCELLS(NAMES, POSITIONS, COLORS, A, COLOR_MARGIN, VIEW_ANGLE,
%           FIG_TITLE, PLOT_TITLE)
%
% Input:
%   names        = the cell names
%   positions    = the cell positions
%   colors       = the cell colors
%   A            = the cell adjacenecy matrix
%   color_margin = the color margin for cell color discrimination
%   view_angle   = the viewing angle for the plot
%   fig_title    = the figure title
%   plot_title   = the plot title

% Are we plotting a line between cells violating the color margin?
A = [];
if length(varargin) > 1
    A = varargin{1};
    color_margin = varargin{2};
end

% Is there a a viewing angle?
view_angle = [0,-90];
if length(varargin) > 2 && ~isempty(varargin{3})
    view_angle = varargin{3};
end

% Is there a figure title?
fig_title = 'Cell Coloring';
if length(varargin) > 3
    fig_title = varargin{4};
end

% Is there a plot title?
plot_title = 'Cells';
if length(varargin) > 4
    plot_title = varargin{5};
end

% Are there too many colors?
color_warning = '';
num_colors = size(colors,2);
if num_colors > 3
    color_warning = ' - WARNING: ONLY SHOWING 3 COLORS!';
end

% Helper function(s).
stretch = @(x,s)([x(1)-s x(2)+s]);

% Initialize the figure.
black = zeros(1,3);
white = ones(1,3);
dark_gray = white/4;
x_pad=3;
y_pad=3;
z_pad=3;
marker_size = 20;
font_size = 8;
font_weight = 'bold';
green_text_thresh = 3/4;
black_edge_thresh = 1/3;

% Plot the cells.
violations = 0;
if ~isempty(fig_title)
    figure('Name', [fig_title color_warning], 'NumberTitle',  'off');
    fig_size = get(0, 'Screensize');
    set(gcf, 'Position', fig_size);
end
for i=1:length(names)
    hold on;
    
    % Determine the color.
    color = zeros(1,3);
    color(1:num_colors) = colors(i,:);
    if length(color) > 3
        color = color(1:3);
    end
    
    % Determine the remaining plot details.
    font_color = white;
    if color(2) > green_text_thresh
        font_color = dark_gray;
    end
    edge_color = black;
    if sum(color) < black_edge_thresh
        edge_color = white;
    end
    
    % Plot the cells.
    plot3(positions(i,1), positions(i,2), positions(i,3), ...
        'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', edge_color, ...
        'MarkerSize', marker_size);
    set(gca, 'Color', 'k', ...
        'Xlim',[min(positions(:,1)) max(positions(:,1))], ...
        'Ylim',[min(positions(:,2)) max(positions(:,2))], ...
        'Zlim',[min(positions(:,3)) max(positions(:,3))]);

    % Show the cell's name.
    text(positions(i,1), positions(i,2), positions(i,3), ...
        names{i}, 'HorizontalAlignment', 'center', ...
        'Color', font_color, 'FontSize', font_size, ...
        'FontWeight', font_weight);
    
    % Plot a line between cells violating the color margin.
    if ~isempty(A)
        for j=1:size(colors,1)
            if A(i,j) > 0 && norm(colors(i,:) - colors(j,:)) < color_margin
                violations = violations + 1;
                x_line = [positions(i,1) positions(j,1)];
                y_line = [positions(i,2) positions(j,2)];
                z_line = [positions(i,3) positions(j,3)];
                plot3(x_line, y_line, z_line, 'w-', 'LineWidth', 2)
            end
        end
    end
end

% Clean up the plot.
axis equal;
axis tight;
set(gca, 'Xlim', stretch(get(gca,'Xlim'), x_pad), ...
    'Ylim', stretch(get(gca,'Ylim'), y_pad), ...
    'Zlim', stretch(get(gca,'Zlim'), z_pad));
view(view_angle);

% Label the plot.
xlabel('Microns');
ylabel('Microns');
zlabel('Microns');
if ~isempty(A)
    plot_title = [plot_title ' - color violations = ' num2str(violations)];
end
title(plot_title);
end
