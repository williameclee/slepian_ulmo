%% GIA2PLMT
% Reads a GIA model and converts it to the spheircal harmonic format.
%
% Syntax
%   Plmt = gia2plmt(model)
%		Returns the GIA geoid change in one year.
%   Plmt = gia2plmt(model, days)
%		Returns the GIA geoid change in the given number of days.
%   Plmt = gia2plmt(model, time)
%		Returns the GIA geoid changes at the given times.
%   Plmt = gia2plmt(model, time, L)
%		Returns the GIA geoid changes truncated to degree L.
%   Plmt = gia2plmt(__, 'Name', Value)
%   [Plmt, PlmtU, PlmtL] = gia2plmt(__)
%		Also returns the upper and lower bounds of the GIA model, if available.
%
% Input arguments
%   model - Name of the GIA model
%		The default model is 'Paulson07'.
%   days - Number of days to calculate the GIA change for
%   time - Vector of dates to calculate the GIA change for
%		The input can be in DATENUM or DATETIME format.
%   L - Maximum degree of the GIA model
%		If empty, the model is not truncated.
%
% Output arguments
%   Plmt - GIA change in spherical harmonic format
%   PlmtU, PlmtL - Upper and lower bounds of the GIA change, if available
%
% See also
%   CORRECT4GIA
%
% Last modified by
%   2024/08/15, williameclee@arizona.edu (@williameclee)

function varargout = gia2plmt(varargin)
    %% Initialisation
    [time, model, L, ~] = parseinputs(varargin{:});

    %% Loading the model
    % Where the model save files are kept
    if strncmp(model, 'Morrow', 6)
        inputFolder = fullfile(getenv('IFILES'), 'GIA', model(1:6));
    else
        inputFolder = fullfile(getenv('IFILES'), 'GIA', model);
    end

    % And the appropriate name
    inputPath = fullfile(inputFolder, sprintf('%s_SD.mat', model));

    % Load this data (saved as lmcosiM)
    load(inputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');
    hasBounds = exist('lmcosiU', 'var') && exist('lmcosiL', 'var');

    %% Some additional processing
    % Truncate the model to the desired degree
    if ~isempty(L)
        lmcosiM = lmcosiM(addmup(L), :);

        if hasBounds
            lmcosiU = lmcosiU(addmup(L), :);
            lmcosiL = lmcosiL(addmup(L), :);
        end

    end

    if isempty(time) || isscalar(time)
        % If time is scalar, interpret it as the day change
        if isempty(time)
            deltaYear = 1;
        else
            deltaYear = time / 365;
        end

        GIAt = lmcosiM;
        GIAt(:, 3:4) = deltaYear * GIAt(:, 3:4);

        if hasBounds
            GIAtU = lmcosiU;
            GIAtU(:, 3:4) = deltaYear * GIAtU(:, 3:4);
            GIAtL = lmcosiL;
            GIAtL(:, 3:4) = deltaYear * GIAtL(:, 3:4);
        end

    else
        % Otherwise, interpret it as a vector of dates (in datenum format)
        % Reference the date string to the first date
        deltaYear = (time - time(1)) / 365;

        GIAt = plm2plmt(lmcosiM, deltaYear);

        if hasBounds
            GIAtU = plm2plmt(lmcosiU, deltaYear);
            GIAtL = plm2plmt(lmcosiL, deltaYear);
        end

    end

    %% Collecting outputs
    if ~hasBounds
        GIAtU = nan;
        GIAtL = nan;

        if nargout > 1
            warning('Upper and lower bounds are not available for this model');
        end

    end

    varargout = {GIAt, GIAtU, GIAtL};

    if nargout > 0
        return
    end

    %% Plotting
    plotgiamap(GIAt, time, deltaYear, model)
end

%% Subfunctions
function varargout = parseinputs(varargin)
    modelD = 'Paulson07';
    p = inputParser;
    addOptional(p, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || isempty(x));
    addOptional(p, 'L', [], ...
        @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'BeQuiet', false, @islogical);

    parse(p, varargin{:});
    time = p.Results.Time;
    model = conddefval(p.Results.Model, modelD);
    L = p.Results.L;
    beQuiet = p.Results.BeQuiet;

    if isdatetime(time)
        time = datenum(time); %#ok<DATNM>
    end

    varargout = {time, model, L, beQuiet};

end

function plmt = plm2plmt(plm, deltaYear)
    plmt = zeros([length(deltaYear), size(plm)]);

    plmt(:, :, 1:2) = repmat(reshape(plm(:, 1:2), [1, length(plm), 2]), [length(deltaYear), 1, 1]);
    plmt(:, :, 3:4) = deltaYear(:) .* reshape(plm(:, 3:4), [1, length(plm), 2]);
    plmt = squeeze(plmt);

end

function plotgiamap(GIAt, time, deltaYear, model)
    % Get the change rate
    if isempty(time) || isscalar(time)
        GIAchange = GIAt;
        GIAchange(:, 3:4) = GIAchange(:, 3:4) / deltaYear;
    else
        GIAchange = squeeze(GIAt(end, :, :) - GIAt(1, :, :));
        GIAchange(:, 3:4) = GIAchange(:, 3:4) / deltaYear(end);
    end

    [GIAmesh, lon, lat] = plm2xyz(GIAchange, "BeQuiet", true);

    % Get the coastlines
    coastLonlat = gshhscoastline('c', 'LonOrigin', 180, "BeQuiet", true);

    [cLim, cStep] = optimalclim(GIAmesh, 'Percentile', 3);
    GIAmesh = max(min(GIAmesh, cLim(2)), cLim(1));

    figure(999)
    set(gcf, "Name", sprintf('GIA change (%s)', upper(mfilename)), "NumberTitle", "off")
    clf

    title(sprintf('Model: %s', model))

    [~, cLevels] = loadcbar(cLim, cStep, "Title", 'GIA change [kg/m^2/yr]', "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, GIAmesh, cLevels, "LineStyle", 'none');
    plotqdm(coastLonlat, 'k');
    hold off
end
