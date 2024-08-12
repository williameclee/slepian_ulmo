%% POLY2XY
%  Converts a polygon to a well-defined curve.
%
% Syntax
%  lonlat = poly2xy(p)
%  lonlat = poly2xy(p, addanchors)
%  [lonlat, lon, lat] = poly2xy(__)
%
% Input arguments
%  p - The polygon to convert
%      The polygon should be a polyshape object.
%  addanchors - Whether to add anchors to the curve
%      The default value is false.
%
% Output arguments
%  lonlat - The curve in the form of a 2-column matrix of longitudes and
%      latitudes.
%  lon, lat - The longitudes and latitudes of the curve
%
% Last modified by
%   2024/08/10, williameclee@arizona.edu (@williameclee)

function [lonlat, lon, lat] = poly2xy(poly, varargin)
    %% Initialisation
    % Parse inputs
    p = inputParser;
    addRequired(p, 'Polygon', @(x) isa(x, 'polyshape'));
    addOptional(p, 'AddAnchors', false, ...
        @(x) islogical(x) || isnumeric(x));
    parse(p, poly, varargin{:});
    poly = p.Results.Polygon;
    anc = p.Results.AddAnchors;

    %% Converting the polygon to a curve
    % Make sure the polygon is closed
    lonlat = closecoastline(poly.Vertices);
    % Make sure the polygon is clockwise
    [lon, lat] = poly2cw(lonlat(:, 1), lonlat(:, 2));
    % Remove duplicate vertices
    lonlat = removeduplicatevertices([lon, lat]);
    % Make sure the minimum longitude is greater than 0
    lonlat(:, 1) = lonlat(:, 1) - floor(min(lonlat(:, 1)) / 360) * 360;

    % Add anchors to the curve
    if anc == 0
        return
    end

    lonlat = addanchors(lonlat);
end
