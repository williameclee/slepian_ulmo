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
%       Structure (1.0) [Data set]. Zenodo. doi: 10.5281/zenodo.5560862.
%
% Last modified by
%   2024/08/16, williameclee@arizona.edu (@williameclee)

function varargout = gia2plmt(varargin)
    %% Initialisation
    [time, model, L, beQuiet] = parseinputs(varargin{:});

    %% Loading the model
    % Load this data (saved as lmcosiM)
    inputPath = finddatafile(model);
    load(inputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');

    if ~beQuiet
        fprintf('%s loaded model %s\n', upper(mfilename), inputPath);
    end

    hasBounds = exist('lmcosiU', 'var') && exist('lmcosiL', 'var');

    %% Some additional processing
    % Truncate the model to the desired degree
    if ~isempty(L)
        lmcosiM = lmcosiM(1:addmup(L), :);

        if hasBounds
            lmcosiU = lmcosiU(1:addmup(L), :);
            lmcosiL = lmcosiL(1:addmup(L), :);
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
            warning('SLEPIAN:gia2plmt:noBoundsToReturn', ...
            'Upper and lower bounds are not available for this model');
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
    modelD = 'Steffen_ice6g_vm5a';

    % Allow skipping the time argument
    if nargin > 0 && ...
            (ischar(varargin{1}) || isstring(varargin{1}) || iscell(varargin{1}))
        varargin(2:end + 1) = varargin(1:end);
        varargin(1) = [];
    end

    p = inputParser;
    addOptional(p, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'L', [], ...
        @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));

    parse(p, varargin{:});
    time = p.Results.Time;
    model = conddefval(p.Results.Model, modelD);
    L = p.Results.L;
    beQuiet = p.Results.BeQuiet;

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

    varargout = {time, model, L, beQuiet};

end

function inputPath = finddatafile(model)

    if exist(model, 'file') == 2
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
    elseif strncmp(model, 'Steffen', 7)
        inputFolder = fullfile(inputFolder, 'SteffenGrids');
    else
        inputFolder = fullfile(inputFolder, model);
    end

    % And the appropriate name
    inputPath = fullfile(inputFolder, sprintf('%s_SD.mat', model));

    if exist(inputPath, 'file') ~= 2
        error('Model %s not found\nIt should be kept at %s', ...
            upper(model), inputPath);
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

function plotgiamap(GIAt, time, deltaYear, model)
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

    [cLim, cStep] = optimalclim(GIAmesh, 'Percentile', 3);
    GIAmesh = max(min(GIAmesh, cLim(2)), cLim(1));

    figure(999)
    % Protect underscore in model name
    model = strrep(model, '_', '\_');
    set(gcf, "NumberTitle", "off", "Name", ...
        sprintf('GIA change in %.1f year(s) (%s)', deltaYear, upper(mfilename)))
    clf

    title(sprintf('Model: %s', model))

    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", 'GIA change [kg/m^2]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, GIAmesh, cLevels, "LineStyle", 'none');
    plotqdm(coastLonlat, 'k');
    hold off
end