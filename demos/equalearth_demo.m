%% EQUALEARTH_DEMO
% This is a demo for EQUALEARTH and EQUALEARTHD.
%
% Last modified by
%   2024/08/14, williameclee@arizona.edu (@williameclee)

function equalearth_demo(varargin)
    %% Initialisation
    if nargin == 0
        funName = '';
    else
        funName = sprintf(' (%s)', upper(char(varargin{1})));
    end

    %% Generating data
    lonOrigin = 200;

    % The shape of the oceans
    lonlat = GeoDomain('oceans').Lonlat("Anchors", true);
    XY = equalearthd(lonlat, lonOrigin);

    % The bounding box
    bbox = [20, -90; 380, -90; 380, 90; 20, 90; 20, -90];
    bbox = addanchors(bbox);
    bboxXY = equalearthd(bbox, lonOrigin);

    %% Plotting
    figName = sprintf('Equal Earth Projection%s', funName);
    figure(999)
    set(gcf, 'Name', figName, 'NumberTitle', 'off')
    clf

    % Plot the original data
    subplot(2, 1, 1)
    title('Original data')

    hold on
    plot(polyshape(bbox), 'FaceColor', 'w')
    plot(polyshape(lonlat))
    hold off

    formatlonticks
    formatlatticks
    axis equal
    xlim([min(bbox(:, 1)), max(bbox(:, 1))])
    ylim([min(bbox(:, 2)), max(bbox(:, 2))])
    set(gca, "Box", 'on', "Layer", 'top', "Color", 'none')

    % Plot the Equal Earth projection
    subplot(2, 1, 2)
    title('Equal Earth projection')

    hold on
    plot(polyshape(bboxXY), 'FaceColor', 'w')
    plot(polyshape(XY))
    hold off

    axis equal
    xlim([min(bboxXY(:, 1)), max(bboxXY(:, 1))])
    ylim([min(bboxXY(:, 2)), max(bboxXY(:, 2))])
    set(gca, "XAxisLocation", 'origin', "YAxisLocation", 'origin', ...
        "TickDir", 'both', "Layer", 'top', "Color", 'none')
end
