%% STERIC2PLMT
% Convert steric sea level products to spherical harmonic (Stokes) coefficients
% Output unit: mm (equivalent to kg/m^2 water)
%
% TODO:
% 	- Make plot
%
% Authored by
%   2025/02/13, williameclee@arizona.edu (@williameclee)

function varargout = steric2plmt(varargin)
    %% Initialisation
    [pcentre, product, L, timelim, beQuiet, outputFormat] = parseinputs(varargin{:});

    %% Loading the model
    % Load this data (saved as lmcosiM)
    [inputSphtPath, inputExists] = finddatafile(pcentre, product, 'spht');

    if inputExists
        load(inputSphtPath, 'stericSpht', 'time');
    else
        [inputLonlattPath, inputExists] = finddatafile(pcentre, product, 'lonlatt');

        if ~inputExists
            error('Model %s not found\nIt should be kept at %s', ...
                upper(product), inputLonlattPath)
        end

        load(inputLonlattPath, 'steric', 'time');

        stericSpht = zeros([addmup(L), 4, length(time)]);

        for iTime = 1:length(time)
            stericSpht(:, :, iTime) = xyz2plm_new(steric(:, :, iTime), L, "BeQuiet", beQuiet);
        end

        save(inputSphtPath, 'stericSpht', 'time')
    end

    %% Collecting outputs
    if ~isempty(timelim)
        isValidTime = time >= timelim(1) & time <= timelim(2);
        time = time(isValidTime);
        stericSpht = stericSpht(:, :, isValidTime);
    end

    switch outputFormat
        case 'timefirst'
            stericSpht = permute(stericSpht, [2, 3, 1]);
        case 'traditional'
            % Do nothing
    end

    varargout = {stericSpht, time};

    if nargout > 0
        return
    end

    %% Plotting
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
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'OutputFormat', 'traditional', @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));

    parse(p, varargin{:});
    pcentre = conddefval(p.Results.Pcenter, pcentreD);
    product = conddefval(p.Results.Product, productD);
    L = conddefval(p.Results.L, LD);
    timelim = p.Results.TimeRange;
    beQuiet = p.Results.BeQuiet;
    outputFormat = p.Results.OutputFormat;

    product = lower(product);

    varargout = {pcentre, product, L, timelim, beQuiet, outputFormat};

end

function [inputPath, inputExists] = finddatafile(pcentre, product, format)

    if ~isempty(getenv('STERIC'))
        inputFolder = getenv('STERIC');
    elseif ~isempty(getenv('IFILES'))
        inputFolder = fullfile(getenv('IFILES'), 'STERIC');
    else
        error('Steric sea level folder not found')
    end

    inputFile = sprintf('%s-%s-%s.mat', pcentre, product, format);
    inputPath = fullfile(inputFolder, inputFile);

    inputExists = exist(inputPath, 'file');

end

function plmt = plm2plmt(plm, deltaYear)
    plmt = zeros([length(deltaYear), size(plm)]);

    plmt(:, :, 1:2) = repmat(reshape( ...
        plm(:, 1:2), [1, length(plm), 2]), [length(deltaYear), 1, 1]);
    plmt(:, :, 3:4) = ...
        deltaYear(:) .* reshape(plm(:, 3:4), [1, length(plm), 2]);
    plmt = squeeze(plmt);

end