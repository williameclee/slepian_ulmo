%% GIA2PLMT
% Reads a GIA model and converts it to the spheircal harmonic format.
%
% Syntax
%   Plmt = gia2plmt(model)
%		Returns the GIA geoid change in one year.
%   Plmt = gia2plmt(days, model)
%		Returns the GIA geoid change in the given number of days.
%   Plmt = gia2plmt(time, model)
%		Returns the GIA geoid changes at the given times.
%   Plmt = gia2plmt(time, model, L)
%		Returns the GIA geoid changes truncated to degree L.
%   Plmt = gia2plmt(__, 'Name', Value)
%   [Plmt, PlmtU, PlmtL] = gia2plmt(__)
%		Also returns the upper and lower bounds of the GIA model, if
%       available.
%
% Input arguments
%   model - Name of the GIA model
%       - 'Steffen_ice6g_vm5a': A model computed by H. Steffen using the
%           ice6g ice model and vm5a viscosity profile. Other models from
%           this dataset are also available and use the original naming
%           scheme. For example, 'Steffen_anu-ice_i72'.
%           This family of models can also be specified as a cell array,
%           e.g. {'Steffen', 'ice6g', 'vm5a'}.
%       - 'LM17.3': A model based on the data from the LM17.3 dataset.
%       - 'Paulson07': A model based on the ICE-5G ice load model of
%           Peltier (2004). Suitable for both Antarctica and Greenland. As
%           corrected by Geruo A and J. Wahr.
%           Please avoid using this model for oceans.
%       - 'Wangetal08': A model based on the older ICE-4G ice model, and
%           viscosity which varies laterally. Suitable for Greenland.
%       - 'IJ05_R2': A model based on the Ivins et al. (2013) ice model.
%           Suitable for Antarctica.
%       - 'IJ05': A model based on the Ivins and James (2005) ice model.
%           Suitable for Antarctica.
%       - 'W12a_v1': A 'best' model from Whitehouse et al. (2012). Suitable
%           only for Antarctica.
%       The input can also be the path to the model file.
%		The default model is 'Steffen_ice6g_vm5a'.
%   days - Number of days to calculate the GIA change for
%   time - Vector of dates to calculate the GIA change for
%		The input can be in DATENUM or DATETIME format.
%   L - Maximum degree of the GIA model
%		If empty, the model is not truncated.
%   OutputField - Field to output
%       - 'massdensity': Surface mass density
%       - 'geoid': Geoid change
%       The default field is 'massdensity'.
%   OutputFormat - Format of the output
%       - 'timefirst': The first dimension is time
%       - 'traditional': The first dimension is degree
%		The default format is 'timefirst'.
%   BeQuiet - Whether to surpress output messages
%		The default option is false.
%
% Output arguments
%   Plmt - GIA change in spherical harmonic format
%   PlmtU, PlmtL - Upper and lower bounds of the GIA change, if available
%
% See also
%   CORRECT4GIA
%
% Data source
%   Paulson, A., S. Zhong, and J. Wahr. (2007). Inference of mantle
%       viscosity from GRACE and relative sea level data, Geophys. J. Int.
%       171, 497–508. doi: 10.1111/j.1365-246X.2007.03556.x
%   Geruo, A., Wahr, J. & Zhong, S. (2013). Computations of the
%       viscoelastic response of a 3-D compressible Earth to surface
%       loading: An application to Glacial Isostatic Adjustment in
%       Antarctica and Canada. Geophys. J. Int. 192, 557-572.
%   Ivins, E. R., T. S. James, J. Wahr, E. J. O. Schrama, F. W. Landerer,
%       and K. M. Simon. (2013). Antarctic contribution to sea level rise
%       observed by GRACE with improved GIA correction, Journal of
%       Geophysical Research: Solid Earth, vol. 118, 3126-3141, doi:
%       10.1002/jgrb.50208
%   Ivins, E. R., T. S. James. (2005). Antarctic glacial isostatic
%       adjustment: a new assessment, Antarctic Science, vol. 17(4),
%       541-553,  doi: 10.1017/S0954102005002968
%   Wang, H., P. Wu, and W. van der Wal. (2018). Using postglacial sea
%       level, crustal velocities and gravity-rate-of-change to constrain
%       the influence of thermal effects on mantle lateral heterogeneities,
%       Journal of Geodynamics 46, 104-117. doi: 10.1016/j.jog.2008.03.003
%   Whitehouse, P. L., Bentley, M. J., Milne, G. A., King, M. A., Thomas,
%       I. D. (2012). A new glacial isostatic adjustment model for
%       Antarctica: calibrated and tested using observations of relative
%       sea-level change and present-day uplift rates. Geophysical Journal
%       International 190, 1464-1482. doi:10.1111/j.1365-246X.2012.05557.x
%   Steffen, H. (2021). Surface Deformations from Glacial Isostatic
%       Adjustment Models with Laterally Homogeneous, Compressible Earth
%       Structure (1.0) [dataset]. Zenodo. doi: 10.5281/zenodo.5560862.
%   Steffen, H., Li, T., Wu, P., Gowan, E. J., Ivins, E., Lecavalier, B.,
%       Tarasov, L., Whitehouse, P. L. (2021). LM17.3 - a global vertical
%       land motion model of glacial isostatic adjustment [dataset].
%       PANGAEA, doi: 10.1594/PANGAEA.932462
%
% Last modified by
%   2025/06/02, williameclee@arizona.edu (@williameclee)

function varargout = gia2plmt(varargin)
    %% Initialisation
    [time, model, L, outputField, outputFormat, beQuiet] = parseinputs(varargin{:});

    %% Loading the model
    % Load this data (saved as lmcosiM)
    inputPath = finddatafile(model);
    warning('off', 'MATLAB:load:variableNotFound');
    load(inputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');
    lmcosiM(1, 3) = 0; % Ensure the mean is zero

    if ~beQuiet
        fprintf('%s loaded model %s\n', upper(mfilename), inputPath);
    end

    hasBounds = exist('lmcosiU', 'var') && exist('lmcosiL', 'var');

    %% Some additional processing
    % Truncate the model to the desired degree
    if ~isempty(L)

        if size(lmcosiM, 1) < addmup(L)

            if ~beQuiet
                warning('SLEPIAN:gia2plmt:truncation', ...
                    'Model %s resolution lower than the requested degree %d', model, L);
            end

            [lmcosiM(1:addmup(L), 2), lmcosiM(1:addmup(L), 1)] = addmon(L);
        else
            lmcosiM = lmcosiM(1:addmup(L), :);
        end

        if hasBounds

            if size(lmcosiU, 1) < addmup(L) || size(lmcosiL, 1) < addmup(L)
                [lmcosiU(1:addmup(L), 2), lmcosiU(1:addmup(L), 1)] = addmon(L);
                [lmcosiL(1:addmup(L), 2), lmcosiL(1:addmup(L), 1)] = addmon(L);
            else
                lmcosiU = lmcosiU(1:addmup(L), :);
                lmcosiL = lmcosiL(1:addmup(L), :);
            end

        end

    end

    % Surface mass density or geoid height
    switch outputField
        case {'massdensity', 'SD'}
            % Do nothing
        case {'geoid', 'POT'}
            % Convert to geoid height
            lmcosiM = plm2pot(lmcosiM, [], [], [], 5);

            if hasBounds
                lmcosiU = plm2pot(lmcosiU, [], [], [], 5);
                lmcosiL = plm2pot(lmcosiL, [], [], [], 5);
            end

    end

    % Time
    if isempty(time) || isscalar(time)
        % If time is scalar, interpret it as the day change
        if isempty(time)
            deltaYear = 1;
        else
            deltaYear = time / days(years(1));
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
        deltaYear = (time - time(1)) / days(years(1));

        GIAt = plm2plmt(lmcosiM, deltaYear);

        if hasBounds
            GIAtU = plm2plmt(lmcosiU, deltaYear);
            GIAtL = plm2plmt(lmcosiL, deltaYear);
        end

    end

    %% Collecting outputs
    switch outputFormat
        case 'timefirst'
            % Do nothing
        case 'traditional'
            GIAt = permute(GIAt, [2, 3, 1]);
    end

    if ~hasBounds
        GIAtU = [];
        GIAtL = [];

        if nargout > 1
            warning('SLEPIAN:gia2plmt:noBoundsToReturn', ...
            'Upper and lower bounds are not available for this model');
        end

    else

        switch outputFormat
            case 'timefirst'
                % Do nothing
            case 'traditional'
                GIAtU = permute(GIAtU, [2, 3, 1]);
                GIAtL = permute(GIAtL, [2, 3, 1]);
        end

    end

    varargout = {GIAt, GIAtU, GIAtL};

    if nargout > 0
        return
    end

    %% Plotting
    plotgiamap(GIAt, time, deltaYear, model, outputField)
end

%% Subfunctions
function varargout = parseinputs(varargin)
    modelD = 'Steffen_ice6g_vm5a';

    % Allow skipping the time argument
    if nargin > 0 && ...
            (ischar(varargin{1}) || isstring(varargin{1}) || iscell(varargin{1}))
        varargin(2:end + 1) = varargin;
        varargin{1} = [];
    end

    p = inputParser;
    addOptional(p, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isempty(x) || idsuration(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'L', [], ...
        @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'massdensity', 'geoid', 'SD', 'POT'})));
    addParameter(p, 'OutputFormat', 'timefirst', ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));

    parse(p, varargin{:});
    time = p.Results.Time;
    model = conddefval(p.Results.Model, modelD);
    L = p.Results.L;
    beQuiet = p.Results.BeQuiet;
    outputField = p.Results.Unit;
    outputFormat = p.Results.OutputFormat;

    if isdatetime(time)
        time = datenum(time); %#ok<DATNM>
    end

    if iscell(model)

        if ~strcmpi(model{1}, 'Steffen')
            error('Unrecognised model name %s', upper(model{1}));
        elseif ~ismember(model{2}, {'anu-ice', 'ice6g', 'ice7g'})
            error('Unrecognised ice model %s', upper(model{2}));
        end

        model = sprintf('%s_%s_%s', model{1}, model{2}, model{3});

    end

    if ismember(lower(model), {'ice6gd', 'ice6g_d', 'ice6g-d', 'ice-6g\_d', 'ice6g\_d'})
        model = 'ICE-6G_D';
    end

    varargout = {time, model, L, outputField, outputFormat, beQuiet};

end

function inputPath = finddatafile(model)

    if isfile(model)
        inputPath = model;
        return
    end

    if ~isempty(getenv('GIA'))
        inputFolder = getenv('GIA');
    elseif ~isempty(getenv('IFILES'))
        inputFolder = fullfile(getenv('IFILES'), 'GIA');
    else
        error('GIA folder not found')
    end

    if strncmp(model, 'Morrow', 6)
        inputFolder = fullfile(inputFolder, model(1:6));
    elseif strncmpi(model, 'Steffen', 7)
        inputFolder = fullfile(inputFolder, 'SteffenGrids');
    elseif strcmp(model, 'LM17.3')
        inputFolder = fullfile(inputFolder, 'LM17.3');

        if exist(inputFolder, 'dir') ~= 7
            mkdir(inputFolder)
        end

    else
        inputFolder = fullfile(inputFolder, model);
    end

    % And the appropriate name
    inputPath = fullfile(inputFolder, sprintf('%s_SD.mat', model));

    if exist(inputPath, 'file') ~= 2

        if strcmp(model, 'LM17.3')

            if exist(fullfile(inputFolder, 'LM17.3_0.5x0.5_geoid_globe.txt'), 'file')
                lm17_Sd(inputFolder, inputPath);
            else
                error('Model %s not found\nPlease download it from %s', ...
                    upper(model), 'https://sites.google.com/view/holgersteffenlm/startseite/data');
            end

        else
            error('Model %s not found\nIt should be kept at %s', ...
                upper(model), inputPath);
        end

    end

end

function plmt = plm2plmt(plm, deltaYear)
    plmt = zeros([length(deltaYear), size(plm)]);

    plmt(:, :, 1:2) = repmat(reshape( ...
        plm(:, 1:2), [1, length(plm), 2]), [length(deltaYear), 1, 1]);
    plmt(:, :, 3:4) = ...
        deltaYear(:) .* reshape(plm(:, 3:4), [1, length(plm), 2]);
    plmt = squeeze(plmt);

end

function plotgiamap(GIAt, time, deltaYear, model, outputField)
    % Get the change rate
    if isempty(deltaYear)
        deltaYear = 1;
    else
        deltaYear = deltaYear(end);
    end

    if isempty(time) || isscalar(time)
        GIAchange = GIAt;
    else
        GIAchange = squeeze(GIAt(end, :, :) - GIAt(1, :, :));
        GIAchange(:, 1:2) = GIAt(1, :, 1:2);
    end

    [GIAmesh, lon, lat] = plm2xyz(GIAchange, "BeQuiet", true);

    % Get the coastlines
    coastLonlat = gshhscoastline('c', 'LonOrigin', 180, "BeQuiet", true);

    [cLim, cStep] = optimalclim(GIAmesh, 'Percentile', 1);
    GIAmesh = max(min(GIAmesh, cLim(2)), cLim(1));

    figure(999)
    % Protect underscore in model name
    model = strrep(model, '_', '\_');
    set(gcf, "NumberTitle", "off", "Name", ...
        sprintf('GIA change in %.1f year(s) (%s)', deltaYear, upper(mfilename)))
    clf

    title(sprintf('Model: %s', model))

    switch outputField
        case {'massdensity', 'SD'}
            cLabel = 'Surface mass density [kg/m^2]';
        case {'geoid', 'POT'}
            cLabel = 'Geoid rate [m/s]';
    end

    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", cLabel, ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, GIAmesh, cLevels, "LineStyle", 'none');
    plotqdm(coastLonlat, 'k');
    hold off
end
