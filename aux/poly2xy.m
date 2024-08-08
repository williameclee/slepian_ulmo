function [lonlat, lon, lat] = poly2xy(poly, varargin)
    p = inputParser;
    addRequired(p, 'Polygon', @(x) isa(x, 'polyshape'));
	addOptional(p, 'Upscale', 0, @(x) isnumeric(x));
	parse(p, poly, varargin{:});
	poly = p.Results.Polygon;
	upscale = p.Results.Upscale;

    lonlat = closecoastline(poly.Vertices);
    [lon, lat] = poly2cw(lonlat(:, 1), lonlat(:, 2));
    lonlat = removeduplicatevertices([lon, lat]);
    % Make sure the minimum longitude is greater than 0
    lonlat(:, 1) = lonlat(:, 1) - floor(min(lonlat(:, 1)) / 360) * 360;

	if upscale == 0
		return
	end

	% Upscale the polygon
	lonlat = addanchors(lonlat);
end
