function varargout = gshhsstruct(varargin)
    %% Initialisation
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    warning('off', 'MATLAB:polyshape:checkAndSimplify');
    p = inputParser;
    addOptional(p, 'DataQuality', 'c', ...
        @(x) ischar(x) && ismember(x, 'cfhil'));
    addOptional(p, 'Upscale', 0, @isnumeric);
    addOptional(p, 'Buffer', 0, @isnumeric);
    addOptional(p, 'MinLandArea', 90 ^ 2, @isnumeric);
    addOptional(p, 'Tolerence', 0, @isnumeric);
    addParameter(p, 'ForceReload', false, @islogical);
    addParameter(p, 'SaveData', true, @islogical);
    addParameter(p, 'Quiet', false, @islogical);
    parse(p, varargin{:});

    dataQuality = p.Results.DataQuality;
    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    minLandArea = p.Results.MinLandArea;
    tol = p.Results.Tolerence;
    forceReload = p.Results.ForceReload;
    saveData = p.Results.SaveData;
    beQuiet = p.Results.Quiet;

    %% Checking if the data already exists
    [matFileName, ~, matFileExists] = gshhsfilename( ...
        'DataQuality', dataQuality, 'Upscale', upscale, ...
        'MinLandArea', minLandArea, 'Tolerence', tol, 'Buffer', buf);

    if matFileExists && ~forceReload
        load(matFileName, 'GshhsCoasts')

        if ~beQuiet
            fprintf('%s loading %s\n', upper(mfilename), matFileName);
        end

        varargout = returndata(nargout, GshhsCoasts);
        return
    end

    %% Loading raw data
    rawFileName = fullfile( ...
        getenv('GSHHS'), ['gshhs_', dataQuality, '.b']);
    GshhsCoasts = gshhs(rawFileName);

    %% Cleaning the data
    GshhsCoasts = GshhsCoasts([GshhsCoasts.Area] >= minLandArea);
    % Remove unnecessary fields
    GshhsCoasts = rmfield(GshhsCoasts, ...
        {'Geometry', 'Container', 'Ancestor', 'FormatVersion', ...
         'Source', 'GSHHS_ID'});
    GshhsCoasts = rmfield(GshhsCoasts, ...
        {'CrossesGreenwich', 'CrossesDateline', 'BoundingBox'});
    % Reformating the data to save space
    numPoints = num2cell(uint8([GshhsCoasts.NumPoints]));
    [GshhsCoasts.NumPoints] = numPoints{:};
    level = num2cell(uint8([GshhsCoasts.Level]));
    [GshhsCoasts.Level] = level{:};
    riverLake = num2cell(logical([GshhsCoasts.RiverLake]));
    [GshhsCoasts.RiverLake] = riverLake{:};
    area = num2cell([GshhsCoasts.Area]); % km^2
    [GshhsCoasts.Area] = area{:};
    areaFull = num2cell([GshhsCoasts.AreaFull]); % km^2
    [GshhsCoasts.AreaFull] = areaFull{:};

    % Remove empty shapes
    GshhsCoasts = GshhsCoasts([GshhsCoasts.NumPoints] > 0);
    % Remove shapes that are not land or Antarctica
    isValidShape = ...
        strcmp({GshhsCoasts.LevelString}, 'land') | ... % is land...
        (strcmp({GshhsCoasts.LevelString}, '') & ...
        [GshhsCoasts.Level] == 5); % or is Antarctica

    GshhsCoasts = GshhsCoasts(isValidShape);
    GshhsCoasts = rmfield(GshhsCoasts, ...
        {'RiverLake', 'LevelString'});

    for iShape = 1:length(GshhsCoasts)
        % Make sure the vertices are valid
        lon = GshhsCoasts(iShape).Lon(:);
        lat = GshhsCoasts(iShape).Lat(:);
        [lon, lat] = poly2cw(lon, lat);

        if iShape == 1 % Make sure Eurasia is intact
            [lat, lon] = flatearthpoly(lat, lon, 90);
        else
            [lat, lon] = flatearthpoly(lat, lon, 0);
        end

        if upscale > 1

            try
                lonlat = bezier([lon(:), lat(:)], upscale);
            catch
                warning(['Upscaling shape ' num2str(iShape), ' failed'])
                lonlat = [lon, lat];
            end

            lon = lonlat(:, 1);
            lat = lonlat(:, 2);

        end

        if tol > 0
            [lat, lon] = reducem(lat, lon, tol);
        end

        if buf > 0
            lonlat = buffercoast([lon, lat], buf);
            lon = lonlat(:, 1);
            lat = lonlat(:, 2);
        end

        GshhsCoasts(iShape).Lon = lon;
        GshhsCoasts(iShape).Lat = lat;
        GshhsCoasts(iShape).NumPoints = ...
            uint16(length(GshhsCoasts(iShape).Lon));
        % Generate polygon
        GshhsCoasts(iShape).Polygon = polyshape(lon, lat);
    end

    GshhsCoasts = GshhsCoasts([GshhsCoasts.Area] >= minLandArea);
    GshhsCoasts = GshhsCoasts([GshhsCoasts.NumPoints] >= 2);

    %% Returning requested data
    varargout = returndata(nargout, GshhsCoasts);

    if ~saveData
        return
    end

    if matFileExists
        save(matFileName, 'GshhsCoasts', '-append')
    else
        save(matFileName, 'GshhsCoasts', '-v7.3')
    end

    if ~beQuiet
        fprintf('%s saving %s\n', upper(mfilename), matFileName)
    end

end

%% Subfunctions
function vOut = returndata(nOut, GshhsCoasts)

    if nOut > 0
        vOut = {GshhsCoasts};
    else
        vOut = {};

        figure(10)
        clf
        hold on

        for iShape = 1:length(GshhsCoasts)
            plot(GshhsCoasts(iShape).Polygon)
            plot(GshhsCoasts(iShape).Lon, GshhsCoasts(iShape).Lat, 'k.-')
        end

        hold off

        axis equal
        axis tight
        grid on
        box on
    end

end
