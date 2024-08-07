function landRegion = ...
        addlandregion(landRegion, lonlim, latlim, varargin)
    p = inputParser;
    addRequired(p, 'LandRegion');
    addRequired(p, 'Lonlim');
    addRequired(p, 'Latlim');
    addOptional(p, 'LongitudeOrigin', 0, @isnumeric);
    parse(p, landRegion, lonlim, latlim, varargin{:});

    landRegion = p.Results.LandRegion;
    lonlim = p.Results.Lonlim;
    latlim = p.Results.Latlim;
    lonOrigin = p.Results.LongitudeOrigin;

    if ~isa(landRegion, 'polyshape')
        isPolygon = false;
        landRegion = polyshape(landRegion);
    else
        isPolygon = true;
    end

    for iLims = 1:size(latlim, 1)
        subtractingRegionLat = [latlim(iLims, 2); latlim(iLims, 2); latlim(iLims, 1); latlim(iLims, 1)];
        subtractingRegionLon = [lonlim(iLims, 1); lonlim(iLims, 2); lonlim(iLims, 2); lonlim(iLims, 1)];
        [subtractingRegionLon, subtractingRegionLat] = ...
            poly2cw(subtractingRegionLon, subtractingRegionLat);
        [subtractingRegionLat, subtractingRegionLon] = ...
            flatearthpoly( ...
            subtractingRegionLat, subtractingRegionLon, ...
            lonOrigin);

        subtractingRegion = polyshape( ...
            subtractingRegionLon, subtractingRegionLat);

        landRegion = union(landRegion, subtractingRegion);
    end

    if ~isPolygon
        landRegion = landRegion.Vertices;
        landRegion = closecoastline(landRegion);
    end

end
