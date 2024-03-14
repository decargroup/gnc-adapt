function plotAxis(Rs, varargin)
    %PLOTAXIS Plot the axis using a triad.
    %
    % Inputs:
    % - Rs: A list of rotation matrices [3 x 3 x N]
    %
    % Options:
    % - axis: axis object
    % - color: colors for a triad [1 x 3]
    % - linewidth: width for each line
    %
    % Editor: Kyungmin John Jung
    % Date: 2024-03-12
    % Lab: DECAR Group
    % Institution: McGill University

    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('axis', [], ...
        @(x) isa(x, 'matlab.graphics.axis.Axes'));
    params.addParameter('color', ['r', 'g', 'b'], ...
        @(x) arrayfun(@(y) iscolor(y), x));
    params.addParameter('alpha', 1, ...
        @(x) isnumeric(x) && x <= 1 && x >= 0);
    params.addParameter('linewidth', 2, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    params.parse(varargin{:});

    axis = params.Results.axis;
    triad_color = params.Results.color;
    alpha = params.Results.alpha;
    linewidth = params.Results.linewidth;

    if isempty(axis)
        gca;
    end

    if length(triad_color) == 1
        triad_color = [triad_color, triad_color, triad_color];
    end

    hold on;

    for i=1:size(Rs, 3)
        C = Rs(:, :, i);
        quiver3(0, 0, 0, C(1, 1), C(1, 2), C(1, 3), ...
            triad_color(1), 'linewidth', linewidth);
        quiver3(0, 0, 0, C(2, 1), C(2, 2), C(2, 3), ...
            triad_color(2), 'linewidth', linewidth);
        quiver3(0, 0, 0, C(3, 1), C(3, 2), C(3, 3), ...
            triad_color(3), 'linewidth', linewidth);
    end

    grid on
    hold off

end

