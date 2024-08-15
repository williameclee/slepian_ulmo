%% OCEANPOLY
% Finds the boundary of the specified ocean(s).
% Note that the boundary does not exclude land regions.
%
% Syntax
%   p = oceanpoly(oceanNames)
%   p = oceanpoly(oceanNames, 'Name', Value)
%   [p, latlim, lonmin] = oceanpoly(__)
%
% Input Arguments
%   oceanNames - The name of the ocean(s) to find the boundary of
%       The names are case-insensitive.
%       The names are as in the LimitsOfOceansAndSeas dataset.
%   latlim - The latitude limits of the ocean boundary
%       The default value is [-90, 90].
%   lonOrigin - The longitude origin of the data
%       The default value is 180.
%   BeQuiet - Suppress the output messages
%       The default value is false.
%
% Output arguments
%   p - The polygon of the ocean boundary
%   latlim, lonmin - The latitude and longitude limits of the ocean 
%       boundary
%
% Data source
%   The ocean boundaries are based on IHO's 'Limits of oceans and seas':
%       International Hydrographic Organization & Sieger, R. (2012).
%       doi: 10.1594/PANGAEA.777975
%
% See also
%   LIMITSOFOCEANSANDSEAS
%
% Last modified by
%   2024/08/12, williameclee@arizona.edu (@williameclee)

function [p, latlim, lonmin] = oceanpoly(oceanNames, varargin)
    %% Initialisation
    % Parse inputs
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

    clear p % Avoid confusion with the output polygon

    %% Finding the ocean boundary
    % Load the ocean boundaries
    LimitsOfOceansAndSeas = limitsofoceansandseas('BeQuiet', beQuiet);
    p = union( ...
        [LimitsOfOceansAndSeas( ...
         ismember({LimitsOfOceansAndSeas.Name}, oceanNames) ...
     ).poly]);

    p = simplify(p);
    [lon, lat] = poly2cw( ...
        p.Vertices(:, 1), p.Vertices(:, 2));

    [lon, lat] = closecoastline(lon, lat);
    [lat, lon] = flatearthpoly(lat, lon, lonOrigin);

    p = polyshape(lon, lat);
    p = removepolarcaps(p, latlim, lonOrigin);
    latlim = ...
        [min(p.Vertices(:, 2)), ...
         max(p.Vertices(:, 2))];
    lonmin = ...
        [min(p.Vertices(:, 1)), ...
         max(p.Vertices(:, 1))];
end
