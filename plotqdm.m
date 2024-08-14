%% PLOTDQM
% Plots a quick and dirty map (on a normal axes) of the given data.
%
% Syntax
%   plotqdm(domain)
%   plotqdm(lon, lat)
%   plotqdm(__, c)
%   plotqdm(__, 'Name', Value)
%   h = plotqdm(__)
%
% Input arguments
%   domain - The geographic domain to plot
%       - A GeoDomain object.
%       - A polyshape object.
%       - An N-by-2 longitude-latitude matrix.
%   lon, lat - Vectors of longitude and latitude.
%   c - Colour of the line.
%   The name-value pair arguments are passed to the plot function.
%
% Output arguments
%   h - Handle to the axes.
%
% Last modified by
%   2024/08/12, williameclee@arizona.edu (@williameclee)

function h = plotqdm(varargin)
    [lon, lat, Inputs] = parselonlatinputs(varargin{:});

    lon = lon - floor(min(lon) / 360) * 360;
    xLim = [min(lon), max(lon)] + [-5, 5];

    if xLim(2) - xLim(1) == 370
        xLim = [0, 360];
    end

    yLim = [min(lat), max(lat)] + [-5, 5];
    yLim = min(max(yLim, -90), 90);

    %% Plotting
    plot(lon, lat, Inputs{:})
    axis equal
    axis tight
    grid on

    xlim(xLim)
    ylim(yLim)

    formatlonticks
    formatlatticks

    set(gca, 'Box', 'on', 'Layer', 'top')

    if nargout == 0
        return
    end

    h = gca;
end