%% EQUALEARTH
% Computes the Equal Earth projection of a set of points.
% The Equal Earth projection is an equal-area pseudocylindrical projection,
% making it suitable for assessing surface mass density-related proceses
% (e.g. the GRACE data we are working with).
% This function only differs from EQUALEARTHD in that it accepts inputs in
% radians instead of degrees.
%
% Syntax
% 	equalearth('demo')
% 	[XY, X, Y] = equalearth(lonlat)
% 	[XY, X, Y] = equalearth(p)
% 	[XY, X, Y] = equalearth(lon, lat)
% 	[XY, X, Y] = equalearth(__, lonOrigin)
% 	equalearth(__)
%
% Input arguments
% 	lonlat - A N-by-2 matrix of longitude-latitude pairs in radians
% 	p - A polyshape object, the vertices are in radians
% 	lon, lat - A vector of longitudes and latitudes in radians
% 	lonOrigin - The origin of the longitude in radians
% 		The default value is 0 (i.e. no shifting)
%
% Output arguments
% 	XY - A N-by-2 matrix of the Equal Earth projection of the input points
% 	X, Y - The X- and Y-coordinates of the Equal Earth projection
%
% Examples
% 	Plot the Equal Earth projection of the world map
% 	>> 	load coastlines
% 	>> 	equalearth(coastlon, coastlat)
%
% Data source
%	The tranformation are based on the article 'The Equal Earth map
%   projection':
%		Šavrič et al. (2019).
%		doi: 0.1080/13658816.2018.1504949
%
% See also
%   EQUALEARTHD
%
% Last modified by
%   2021/08/14, williameclee@arizona.edu (@williameclee)

function varargout = equalearth(varargin)
    %% Initialisation
    % Add path to the auxiliary functions
    addpath(fullfile(fileparts(mfilename('fullpath')), 'demos'));

    % Demos
    if strcmpi(varargin{1}, 'demo')
        equalearth_demo(mfilename)
        return
    end

    % Parse inputs
    [lon, lat, lonOrigin] = parseinputs(varargin{:});
    lon = lon - lonOrigin;

    %% Coordinate transformation
    theta = asin(sqrt(3) / 2 * sin(lat));
    lambda = lon;
    A = [1.340264, -0.081106, 0.000893, 0.003796];

    X = (2 * sqrt(3) * lambda .* cos(theta));
    X = X ./ (3 * (9 * A(4) * theta .^ 8 + 7 * A(3) * theta .^ 6 + ...
        3 * A(2) * theta .^ 2 + A(1)));

    Y = A(4) * theta .^ 9 + A(3) * theta .^ 7 + A(2) * theta .^ 3 + ...
        A(1) * theta;

    %% Collecting output
    X = X(:);
    Y = Y(:);
    XY = [X, Y];

    if nargout == 0
        figName = sprintf('Equal Earth Projection (%s)', upper(mfilename));
        figure(999)
        set(gcf, 'Name', figName, 'NumberTitle', 'off')

        plot(X, Y, 'k', 'LineWidth', 0.5)

        axis equal
        xlim([min(X) - 0.1, max(X) + 0.1])
        ylim([min(Y) - 0.1, max(Y) + 0.1])

        return
    end

    if isa(varargin{1}, 'polyshape')
        varargout = {polyshape(XY), X, Y};
    else
        varargout = {XY, X, Y};
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Find longitude and latitude
    [lon, lat, varargin] = parselonlatinputs(varargin{:});

    % Find the longitude origin
    lonOriginD = 180 * double(min(lon) >= 0);
    p = inputParser;
    addOptional(p, 'LonOrigin', lonOriginD, ...
        @(x) isnumeric(x) || ischar(x) || isstring(x));
    parse(p, varargin{:});
    lonOrigin = p.Results.LonOrigin;

    if ischar(lonOrigin) || isstring(lonOrigin)

        if ismember(lonOrigin, {'center', 'centre', 'c'})
            lonOrigin = (max(lon) + min(lon)) / 2;
        else
            error('Invalid lonOrigin argument %s', lonOrigin)
        end

    end

    varargout = {lon, lat, lonOrigin};
end
