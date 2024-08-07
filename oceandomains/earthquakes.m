function XY = earthquakes(varargin)
    p = inputParser;
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 10, ...
        @(x) isnumeric(x));
    addOptional(p, 'Epiccenter', [], @(x) (isnumeric(x) && length(x) == 2) || isempty(x));
    parse(p, varargin{:});

    upscale = p.Results.Upscale;
    eqMaskSize = p.Results.Buffer;
    eqLonlat = p.Results.Epiccenter;

    if isempty(upscale)
        upscale = 0;
    end

    if isempty(eqLonlat)
        eqLonlat = ...
            [95.854, 3.316; ... % 2004 Indian Ocean earthquake
             142.269, 38.322]; % 2011 Tōhoku earthquake
    end

    eqMask = cell(1, size(eqLonlat, 1));

    for i = 1:size(eqLonlat, 1)
        [maskLat, maskLon] = bufferm(eqLonlat(i, 2), eqLonlat(i, 1), eqMaskSize, "outPlusInterior");

        if upscale ~= 0 && upscale ~= 1
            maskLonlat = bezier([maskLon, maskLat], upscale);
        else
            maskLonlat = [maskLon, maskLat];
        end

        eqMask{i} = polyshape(maskLonlat);
    end

    eqMask = union(eqMask{:});
    XY = eqMask.Vertices;

    XY = closecoastline(XY);

    if nargout == 0
        load coastlines %#ok<LOAD>

        figure(10)
        clf

        hold on
        plot(coastlon, coastlat);
        plot(XY(:, 1), XY(:, 2));
        hold off

        axis equal
        grid on

        xlim([min(XY(:, 1)) - 5, max(XY(:, 1)) + 5])
        ylim([min(XY(:, 2)) - 5, max(XY(:, 2)) + 5])
    end

end
