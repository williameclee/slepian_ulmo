%% AOD1B2PLMT
% Fetch the AOD1B modelled GRACE geoid (or surface mass density) coefficient time series
%
% Syntax
%   [dates, dataSpht, stdSpht] = aod1b2plmt(Pcenter, Rlevel, product, L, units)
%
% Inputs
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data center is 'CSR'.
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06'.
%       The default release level is 'RL06'.
%   Product - Type of product
%       - 'GAC'
%       - 'GAD'
%       The default product is 'GAC'.
%   L - The bandwidth of the output product
%       The default L is 60.
%   Unit - Unit of the output
%       - 'POT': Geopotential field.
%       - 'SD': Surface mass density.
%       The default field is 'SD'.
%   TimeRange - Time range of the output
%       When specified, the output will be limited to the specified time range.
%       The default time range is [] (all available data).
%   ForceNew - Whether to force new generation of a save file
%       The default option is false.
%   BeQuiet - Whether to suppress output
%       The default option is false.
%   OutputFormat - Output format
%       - 'timefirst': The output is in time-first format.
%       - 'traditional': The output is in traditional format.
%       See below for details.
%       The default output format is 'timefirst', which is consistent with the new implementations of, e.g. GRACE2PLMT_NEW; the original format in SLEPIAN_ALPHA is 'traditional'.
%
% Outputs
%   dates - Time stamps of data
%       The format is MATLAB's DATETIME.
%   spht - Desired field
%       - 'timefirst': The output is in [nMonths, addmup(L), 4] format.
%       - 'traditional': The output is in [addmup(L), 4, nMonths] format.
%           This is the same format as the original GRACE2PLMT.
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 6]
%                    these could also be in surface mass density
%    cal_errors     calibrated errors [nmonths x addmup(Ldata) x 4]
%
%   See also
%       GRACE2PLMT
%
% Last modified by
%   2025/03/18, williameclee@arizona.edu (@williameclee)
%   2014/02/27, charig@princeton.edu
%   2011/05/17, fjsimons@alum.mit.edu

function varargout = aod1b2plmt(varargin)
    %% Initialisation
    % Parse inputs
    [Pcenter, Rlevel, product, Loutput, units, ...
         timelim, forcenew, beQuiet, outputFormat] = ...
        parseinputs(varargin);

    % If this file already exists, load it.  Otherwise, or if we force it, make
    % a new one (e.g. you added extra months to the database).
    [outputPath, ~, outputExists] = getoutputpath(Rlevel, Pcenter, product, units);

    if outputExists && ~forcenew
        load(outputPath, 'potcoffs', 'dates', 'cal_errors');

        [dates, potcoffs, cal_errors] = ...
            formatoutput(dates, potcoffs, cal_errors, Loutput, timelim, outputFormat);

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), outputPath)
        end

        varargout = {dates, potcoffs, cal_errors};
        return
    end

    %% Computing data
    if ~isempty(getenv('ORIGINALGRACEDATA'))
        inputFolder = fullfile(getenv('ORIGINALGRACEDATA'), ...
            Rlevel, Pcenter);
    elseif ~isempty(getenv('GRACEDATA'))
        inputFolder = fullfile(getenv('GRACEDATA'), ...
            'raw', Rlevel, Pcenter);
    else
        inputFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'raw', Rlevel, Pcenter);
    end

    switch Pcenter
        case 'GFZ'
            error('GFZ not implemented yet')
        case 'CSR'

            switch Rlevel
                case 'RL06'
                    inputFileList = ...
                        ls2cell(fullfile(inputFolder, sprintf('%s-2_*', product)));
                    Ldata = 180;
            end

        case 'JPL'
            error('JPL not implemented yet')
    end

    % Preallocate
    nMonths = length(inputFileList);
    dates = NaT([1, nMonths]);
    potcoffs = nan([addmup(Ldata), 4, nMonths]);
    cal_errors = nan([addmup(Ldata), 4, nMonths]);

    % Loop over the months
    for iMonth = 1:nMonths
        % load geopotential coefficients
        inputPath = fullfile(inputFolder, inputFileList{iMonth});
        [dates(iMonth), potcoffs(:, :, iMonth), cal_errors(:, :, iMonth)] = ...
            aod1b2plm(inputPath, Pcenter, "DateFormat", 'datetime');

        fprintf('%s processed %s (month %3i) \n', upper(mfilename), ...
            datetime(dates(iMonth), "Format", "uuuu/MM/dd HH:mm"), iMonth);
    end

    % Convert to surface mass density
    if strcmp(units, 'SD')
        potcoffs = plm2pot(potcoffs, [], [], [], 4);
        cal_errors = plm2pot(cal_errors, [], [], [], 4);
    end

    % Save
    save(outputPath, 'potcoffs', 'cal_errors', 'dates');

    if ~beQuiet
        fprintf('%s saved %s\n', upper(mfilename), outputPath)
    end

    % Collect output
    [dates, potcoffs, cal_errors] = formatoutput(dates, potcoffs, cal_errors, Loutput, timelim, outputFormat);

    varargout = {dates, potcoffs, cal_errors};
end

%% Subfunctions
function varargout = parseinputs(varargin)
    ip = inputParser;
    addOptional(ip, 'Pcenter', 'CSR', ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', 'RL06', ...
        @(x) ischar(validatestring(x, {'RL04', 'RL05', 'RL06'})));
    addOptional(ip, 'Product', 'GAC', ...
        @(x) ischar(validatestring(x, {'GAC', 'GAD'})));
    addOptional(ip, 'L', 60, ...
        @(x) isnumeric(x) && x > 0);
    addOptional(ip, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'POT', 'SD'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) ((isnumeric(x) || isdatetime(x)) && length(x) == 2) || ...
        isempty(x));
    addOptional(ip, 'ForceNew', false, ...
        @(x) isnumeric(x) || islogical(x));
    addOptional(ip, 'BeQuiet', false, ...
        @(x) isnumeric(x) || islogical(x));
    addParameter(ip, 'OutputFormat', 'timefirst', ... % to be consistent with grace2plmt_new
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));

    varargin = varargin{:};

    if iscell(varargin{1})
        varargin = [varargin{1}{:}, varargin(2:end)];
    end

    parse(ip, varargin{:});

    Pcenter = ip.Results.Pcenter;
    Rlevel = ip.Results.Rlevel;
    product = ip.Results.Product;
    Loutput = ip.Results.L;
    unit = ip.Results.Unit;
    timelim = ip.Results.TimeRange;
    forcenew = logical(ip.Results.ForceNew);
    beQuiet = logical(ip.Results.BeQuiet);
    outputFormat = ip.Results.OutputFormat;

    if isnumeric(timelim) && ~isempty(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    varargout = {Pcenter, Rlevel, product, Loutput, unit, timelim, forcenew, beQuiet, outputFormat};
end

function [outputPath, outputFolder, outputExists] = getoutputpath(Rlevel, Pcenter, product, units)

    if ~isempty(getenv('GRACEDATA'))
        outputFolder = fullfile(getenv('GRACEDATA'));
    else
        outputFolder = fullfile(getenv('IFILES'), 'GRACE');
    end

    % And the name of that save file
    if strcmp(units, 'SD')
        outputPath = fullfile(outputFolder, ...
            sprintf('%s_%s_%s_SD.mat', Pcenter, Rlevel, product));
    else
        outputPath = fullfile(outputFolder, ...
            sprintf('%s_%s_%s.mat', Pcenter, Rlevel, product));
    end

    outputExists = exist(outputPath, 'file');

end

function [dates, potcoffs, cal_errors] = ...
        formatoutput(dates, potcoffs, cal_errors, Loutput, timelim, outputFormat)

    if isnumeric(dates)
        dates = datetime(dates, "ConvertFrom", 'datenum');
    end

    if ~isempty(timelim)
        isValidTime = dates >= timelim(1) & dates <= timelim(2);
        dates = dates(isValidTime);
        potcoffs = potcoffs(:, :, isValidTime);
        cal_errors = cal_errors(:, :, isValidTime);
    end

    potcoffs(3, 1, potcoffs(3, 1, :) == 1839) = 1; % Some qeird bug?
    potcoffs(4, 1, potcoffs(4, 1, :) == 3678) = 2; % Some qeird bug?

    if ~isempty(Loutput)
        potcoffs = potcoffs(1:addmup(Loutput), :, :);
        cal_errors = cal_errors(1:addmup(Loutput), :, :);
    end

    if strcmp(outputFormat, 'timefirst') && size(potcoffs, 1) ~= length(dates)
        potcoffs = permute(potcoffs, [3, 1, 2]);
        cal_errors = permute(cal_errors, [3, 1, 2]);
    elseif strcmp(outputFormat, 'traditional') && size(potcoffs, 3) ~= length(dates)
        potcoffs = permute(potcoffs, [2, 3, 1]);
        cal_errors = permute(cal_errors, [2, 3, 1]);
    end

end
