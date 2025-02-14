%% STERIC2PLMT
% Convert steric sea level products to spherical harmonic (Stokes) coefficients
% Output unit: mm (equivalent to kg/m^2 water)
%
% TODO:
% 	- Make plot
%
% Data Source:
%   - NOAA Global Ocean Heat and Salt Content: Seasonal, Yearly, and Pentadal Fields
%       https://www.ncei.noaa.gov/access/global-ocean-heat-content/
%
% Authored by
%   2025/02/13, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/02/14, williameclee@arizona.edu (@williameclee)

function varargout = steric2plmt(varargin)
    %% Initialisation
    [pcentre, product, L, timelim, forceNew, saveData, beQuiet, outputFormat] ...
        = parseinputs(varargin{:});

    %% Loading the model
    % Load this data (saved as lmcosiM)
    [inputSphtPath, inputExists] = findoutputfile(pcentre, product, L);

    if inputExists && ~forceNew
        load(inputSphtPath, 'stericSpht', 'time', 'hasData');
    else
        [inputLonlattPath, inputExists] = findinputfile(pcentre, product);

        if ~inputExists
            error('Model %s not found\nIt should be kept at %s', ...
                upper(product), inputLonlattPath)
        end

        load(inputLonlattPath, 'steric', 'time');
        hasData = true(size(time));

        if size(steric, 1) > size(steric, 2)
            steric = permute(steric, [2, 1, 3]);
        end

        stericSpht = zeros([addmup(L), 4, length(time)]);

        for iTime = 1:length(time)

            if any(isnan(steric(:, :, iTime)), 'all')
                hasData(iTime) = false;
                continue
            end

            stericSpht(:, :, iTime) = ...
                xyz2plm_new(steric(:, :, iTime), L, "BeQuiet", beQuiet >= 1);
        end

        if saveData
            save(inputSphtPath, 'stericSpht', 'time', 'hasData')
        end

    end

    %% Collecting outputs
    if ~isempty(timelim)
        isValidTime = time >= timelim(1) & time <= timelim(2);
        time = time(isValidTime);
        stericSpht = stericSpht(:, :, isValidTime);
        hasData = hasData(isValidTime);
    end

    switch outputFormat
        case 'timefirst'
            stericSpht = permute(stericSpht, [2, 3, 1]);
        case 'traditional'
            % Do nothing
    end

    varargout = {stericSpht, time, hasData};

    if nargout > 0
        return
    end

    %% Plotting
    [stericLonlat, lon, lat] = plm2xyz(stericSpht(:, :, end), "BeQuiet", beQuiet);
    figure(999)
    clf

    contourf(lon, lat, stericLonlat)
    colorbar
    axis equal
    xlim([0, 360])
    ylim([-90, 90])
end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Default values
    pcentreD = 'NOAA';
    productD = '2000total';
    LD = 60;
    p = inputParser;
    addOptional(p, 'Pcenter', pcentreD, ...
        @(x) ischar(validatestring(x, {'NOAA'})));
    addOptional(p, 'Product', productD, ...
        @(x) ischar(validatestring(x, ...
        {'2000total', '2000thermosteric', '2000halosteric', '700total', '700thermosteric', '700halosteric'})));
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'TimeRange', [], ...
        @(x) (isnumeric(x) || isdatetime(x) && length(x) == 2) || isempty(x));
    addParameter(p, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveData', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeQuiet', 0.5, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'OutputFormat', 'traditional', @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));

    parse(p, varargin{:});
    pcentre = conddefval(p.Results.Pcenter, pcentreD);
    product = conddefval(p.Results.Product, productD);
    L = conddefval(p.Results.L, LD);
    timelim = p.Results.TimeRange;
    forceNew = uint8(p.Results.ForceNew * 2);
    saveData = logical(p.Results.SaveData);
    beQuiet = logical(p.Results.BeQuiet);
    outputFormat = p.Results.OutputFormat;

    product = lower(product);

    varargout = {pcentre, product, L, timelim, forceNew, saveData, beQuiet, outputFormat};

end

function [inputPath, inputExists] = findinputfile(pcentre, product)

    if ~isempty(getenv('STERIC'))
        inputFolder = getenv('STERIC');
    elseif ~isempty(getenv('IFILES'))
        inputFolder = fullfile(getenv('IFILES'), 'STERIC');
    else
        error('Steric sea level folder not found')
    end

    inputFile = sprintf('%s-%s-lonlatt.mat', pcentre, product);
    inputPath = fullfile(inputFolder, inputFile);

    inputExists = exist(inputPath, 'file');

end

function [inputPath, inputExists] = findoutputfile(pcentre, product, L)

    if ~isempty(getenv('STERIC'))
        inputFolder = fullfile(getenv('STERIC'), 'SPH');
    elseif ~isempty(getenv('IFILES'))
        inputFolder = fullfile(getenv('IFILES'), 'STERIC', 'SPH');
    else
        error('Steric sea level folder not found')
    end

    if ~exist(inputFolder, 'dir')
        mkdir(inputFolder)
    end

    inputFile = sprintf('%s-%s-L%d.mat', pcentre, product, L);
    inputPath = fullfile(inputFolder, inputFile);

    inputExists = exist(inputPath, 'file');

end
