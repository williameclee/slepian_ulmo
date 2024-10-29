function varargout = gshhscoastline(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    % Parse inputs
    [dataQuality, latlim, lonlim, minLandArea, ...
         upscale, buf, ~, lonOrigin, ...
         gshhsFileName, gshhsFileExists, ...
         forcenew, saveData, beQuiet] = parseinputs(varargin);

    %% Reteiving the original data
    if gshhsFileExists && ~forcenew
        load(gshhsFileName, ...
            'GshhsCoasts', 'gshhsCoastXY', 'gshhsCoastPoly')

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), gshhsFileName);
        end

    else
        GshhsCoasts = gshhsstruct('DataQuality', dataQuality, ...
            'Upscale', upscale, 'Buffer', buf, ...
            'Quiet', true, 'ForceReload', forcenew);
    end

    % Check if the data is complete
    needsUpdate = false;

    if ~exist('gshhsCoastPoly', 'var')
        needsUpdate = true;

        if ~beQuiet
            fprintf('%s generating the coasline polygon, this may take a while...\n', ...
                upper(mfilename))
        end

        segmentsToInclude = 1:length(GshhsCoasts);
        segmentsToInclude = ...
            segmentsToInclude(all(segmentsToInclude ~= [57, 77, 78]'));
        gshhsCoastPoly = union([GshhsCoasts(segmentsToInclude).Polygon]);
    end

    if ~exist('gshhsCoastXY', 'var')
        gshhsCoastXY = poly2xy(gshhsCoastPoly);
    end

    % Save the data if requested
    if saveData && needsUpdate

        if gshhsFileExists
            save(gshhsFileName, ...
                'gshhsCoastXY', 'gshhsCoastPoly', 'GshhsCoasts', ...
            '-append')

            if ~beQuiet
                fprintf('%s updated %s\n', upper(mfilename), gshhsFileName)
            end

        else
            save(gshhsFileName, ...
                'gshhsCoastXY', 'gshhsCoastPoly', 'GshhsCoasts')

            if ~beQuiet
                fprintf('%s saved %s\n', upper(mfilename), gshhsFileName)
            end

        end

    end

    %% Cropping the data to the desired limits
    % Leave only the land with a minimum area
    GshhsCoasts = GshhsCoasts([GshhsCoasts.Area] >= minLandArea);

    % Crop the data to the desired limits
    gshhsCoastPoly = croptolims(gshhsCoastPoly, ...
        latlim, lonlim, lonOrigin);

    % Make sure the coastline is closed
    gshhsCoastXY = closecoastline(gshhsCoastPoly.Vertices);

    %% Returning requested data
    varargout = {gshhsCoastXY, gshhsCoastPoly, GshhsCoasts};

    if nargout > 0
        return
    end

    figure(10)
    set(gcf, "Name", sprintf('Coastline (%s)', upper(mfilename)), ...
        "NumberTitle", "off")
    clf

    plotqdm(gshhsCoastXY, 'k.-')

end

%% Subfunctions
function varargout = parseinputs(inputArguments)
    %% Defining the input parser
    p = inputParser;
    addOptional(p, 'dataQuality', 'c', ...
        @(x) ischar(x) && ismember(x, 'cfhil'));
    addOptional(p, 'LatLim', [-90, 90], ...
        @(x) isnumeric(x) && length(x) == 2);
    addOptional(p, 'LonLim', [0, 360], ...
        @(x) isnumeric(x) && length(x) == 2);
    addOptional(p, 'MinLandArea', 30, @isnumeric);
    addOptional(p, 'Upscale', 0, @isnumeric);
    addOptional(p, 'Buffer', 0, @isnumeric);
    addOptional(p, 'Tolerence', 0.2, @isnumeric);
    addParameter(p, 'LonOrigin', 180, @isnumeric);
    addParameter(p, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveData', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    parse(p, inputArguments{:});

    %% Assigning the inputs
    dataQuality = p.Results.dataQuality;
    latlim = p.Results.LatLim;
    lonlim = p.Results.LonLim;
    minLandArea = p.Results.MinLandArea;
    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    tol = p.Results.Tolerence;
    lonOrigin = p.Results.LonOrigin;
    forcenew = logical(p.Results.ForceNew);
    saveData = logical(p.Results.SaveData);
    beQuiet = logical(p.Results.BeQuiet);

    %% Additional parameters
    [gshhsFileName, ~, gshhsFileExists] = gshhsfilename( ...
        'DataQuality', dataQuality, ...
        'Upscale', upscale, 'Buffer', buf, ...
        'MinLandArea', minLandArea, 'Tolerence', tol);

    %% Returning the inputs
    varargout = ...
        {dataQuality, latlim, lonlim, minLandArea, ...
         upscale, buf, tol, lonOrigin, ...
         gshhsFileName, gshhsFileExists, ...
         forcenew, saveData, beQuiet};
end

function p = croptolims(p, latlim, lonlim, lonOrigin)
    %% Shifting
    XY = poly2xy(p);
    [Y, X] = ...
        flatearthpoly(XY(:, 2), XY(:, 1), lonOrigin);
    X = X - floor(min(X(:))/360)*360;
    p = polyshape([X, Y]);

    %% Cropping
    bbox = polyshape( ...
        [lonlim(1), lonlim(2), lonlim(2), lonlim(1), lonlim(1)], ...
        [latlim(1), latlim(1), latlim(2), latlim(2), latlim(1)]);
    p = intersect(p, bbox);
end
