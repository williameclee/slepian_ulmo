function varargout = gshhscoastline(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    % Parse inputs
    [dataQuality, latlim, lonlim, minLandArea, ...
         upscale, buf, ~, lonOrigin, ...
         gshhsFileName, gshhsFileExists, ...
         forceReload, saveData, beQuiet] = parseinputs(varargin);

    %% Reteiving the original data
    if gshhsFileExists && ~forceReload
        load(gshhsFileName, ...
            'GshhsCoasts', 'gshhsCoastXY', 'gshhsCoastPoly')

        if ~beQuiet
            fprintf('%s loading %s\n', upper(mfilename), gshhsFileName);
        end

    else
        GshhsCoasts = gshhsstruct('DataQuality', dataQuality, ...
            'Upscale', upscale, 'Buffer', buf, ...
            'Quiet', true, 'ForceReload', forceReload);
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
        segmentsToInclude = segmentsToInclude(all(segmentsToInclude ~= [57, 77, 78]'));
        gshhsCoastPoly = union([GshhsCoasts(segmentsToInclude).Polygon]);
    end

    if ~exist('gshhsCoastXY', 'var')
        gshhsCoastXY = closecoastline(gshhsCoastPoly.Vertices);
    end

    % Save the data if requested
    if saveData && needsUpdate

        if gshhsFileExists
            save(gshhsFileName, ...
                'gshhsCoastXY', 'gshhsCoastPoly', 'GshhsCoasts', ...
            '-append')
        else
            save(gshhsFileName, ...
                'gshhsCoastXY', 'gshhsCoastPoly', 'GshhsCoasts')
        end

        if ~beQuiet
            fprintf('%s saving or updating %s\n', upper(mfilename), gshhsFileName)
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
    if nargout > 0
        varargout = {gshhsCoastXY, gshhsCoastPoly, GshhsCoasts};
    else
        figure(10)
        clf

        plot(gshhsCoastXY(:, 1), gshhsCoastXY(:, 2), ...
        'k.-')

        axis equal
        axis tight
        grid on
        box on
    end

end

%% Subfunctions
function varargout = parseinputs(inputArguments)
    %% Defining the input parser
    p = inputParser;
    addOptional(p, 'dataQuality', 'c', ...
        @(x) ischar(x) && ismember(x, 'cfhil'));
    addOptional(p, 'LatLim', [-90, 90], ...
        @(x) isnumeric(x) && length(x) == 2);
    addOptional(p, 'LonLim', [-180, 180], ...
        @(x) isnumeric(x) && length(x) == 2);
    addOptional(p, 'MinLandArea', 30, @isnumeric);
    addOptional(p, 'Upscale', 0, @isnumeric);
    addOptional(p, 'Buffer', 0, @isnumeric);
    addOptional(p, 'Tolerence', 0.2, @isnumeric);
    addParameter(p, 'LongitudeOrigin', 0, @isnumeric);
    addParameter(p, 'ForceReload', false, @islogical);
    addParameter(p, 'SaveData', true, @islogical);
    addParameter(p, 'Quiet', false, @islogical);
    parse(p, inputArguments{:});

    %% Assigning the inputs
    dataQuality = p.Results.dataQuality;
    latlim = p.Results.LatLim;
    lonlim = p.Results.LonLim;
    minLandArea = p.Results.MinLandArea;
    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    tol = p.Results.Tolerence;
    lonOrigin = p.Results.LongitudeOrigin;
    forceReload = p.Results.ForceReload;
    saveData = p.Results.SaveData;
    beQuiet = p.Results.Quiet;

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
         forceReload, saveData, beQuiet};
end

function gshhsCoastPoly = ...
        croptolims(gshhsCoastPoly, latlim, lonlim, lonOrigin)
    %% Making sure the coastline is well defined
    gshhsCoastXY = gshhsCoastPoly.Vertices;
    gshhsCoastXY = closecoastline(gshhsCoastXY);
    % Specify the longitude origin
    [gshhsCoastY, gshhsCoastX] = ...
        flatearthpoly(gshhsCoastXY(:, 2), gshhsCoastXY(:, 1), lonOrigin);
    gshhsCoastPoly = polyshape(gshhsCoastX, gshhsCoastY);

    %% Cropping the data to the desired limits
    boundingBox = polyshape( ...
        [lonlim(1), lonlim(2), lonlim(2), lonlim(1), lonlim(1)], ...
        [latlim(1), latlim(1), latlim(2), latlim(2), latlim(1)]);
    gshhsCoastPoly = intersect(gshhsCoastPoly, boundingBox);
end
