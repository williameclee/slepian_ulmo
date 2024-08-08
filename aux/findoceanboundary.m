function [oceanLimit, latlim, lonmin] = ...
        findoceanboundary(oceanNames, varargin)
    p = inputParser;
    addRequired(p, 'Oceans', @(x) iscell(x) || ischar(x) || isstring(x));
    addOptional(p, 'Latlim', [-90, 90], @isnumeric);
    addOptional(p, 'LonOrigin', 180, @isnumeric);
    addParameter(p, 'BeQuiet', false, ...
        @(x) islogical(x) || isnumeric(x));
    parse(p, oceanNames, varargin{:});
    oceanNames = p.Results.Oceans;
    latlim = p.Results.Latlim;
    lonOrigin = p.Results.LonOrigin;
    beQuiet = logical(p.Results.BeQuiet);

    LimitsOfOceansAndSeas = limitsofoceansandseas('BeQuiet', beQuiet);
    oceanLimit = union( ...
        [LimitsOfOceansAndSeas( ...
         ismember({LimitsOfOceansAndSeas.Name}, oceanNames) ...
     ).poly]);

    oceanLimit = simplify(oceanLimit);
    [lon, lat] = poly2cw( ...
        oceanLimit.Vertices(:, 1), oceanLimit.Vertices(:, 2));

    [lon, lat] = closecoastline(lon, lat);
    [lat, lon] = flatearthpoly( ...
        lat, lon, lonOrigin);

    oceanLimit = polyshape(lon, lat);
    oceanLimit = removepolarcaps(oceanLimit, latlim, lonOrigin);
    latlim = ...
        [min(oceanLimit.Vertices(:, 2)), ...
         max(oceanLimit.Vertices(:, 2))];
    lonmin = ...
        [min(oceanLimit.Vertices(:, 1)), ...
         max(oceanLimit.Vertices(:, 1))];
end
