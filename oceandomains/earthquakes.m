function XY = earthquakes(varargin)
    p = inputParser;
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 10, ...
        @(x) isnumeric(x));
    addOptional(p, 'Epicenter', [], ...
        @(x) (isnumeric(x) && length(x) == 2) || isempty(x));
    parse(p, varargin{:});

    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    epicentre = p.Results.Epicenter;

    if isempty(upscale) || upscale == 1
        upscale = 0;
    end

    if isempty(epicentre)
        epicentre = ...
            [95.854, 3.316; ... % 2004 Indian Ocean earthquake
             142.269, 38.322]; % 2011 Tōhoku earthquake
    end

    p = cell(1, size(epicentre, 1));

    for i = 1:size(epicentre, 1)
        [maskLat, maskLon] = bufferm(epicentre(i, 2), epicentre(i, 1), ...
            buf, "outPlusInterior");
        p{i} = polyshape(maskLon, maskLat);
    end

    p = union(p{:});
    XY = p.Vertices;

    XY = closecoastline(XY);

    if upscale > 3
        XY = bezier(XY, upscale);
    else
        XY = bezier(XY, 3);
    end

    if nargout > 0
        return
    end

    %% Plotting the boundary
    load coastlines %#ok<LOAD>

    figName = sprintf('Earthquakes (%s)', upper(mfilename));

    figure(999)
    set(gcf, 'Name', figName, 'NumberTitle', 'off')
    clf

    dirtymap(XY)
    hold on
    plot(coastlon, coastlat, 'k')
    plot(XY(:, 1), XY(:, 2), 'b.-')
    hold off
end
