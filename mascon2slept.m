%% MASCON2SLEPT
% Computes Slepian expansions of the CSR RL06 GRACE mascons
%
% Syntax
%   [slept, dates] = mascon2slept(dataType, domain, L)
%   [slept, dates] = mascon2slept(__, timerange)
%   [__, G, V, N, timeNum] = mascon2slept(__)
%
% Input arguments
%   dataType - The type of data to use
%       It can be either 'mascon', 'masconC20', 'masconC30', 'slrC20', 'slrC30', or 'deg123'.
%       The default value is 'mascon' (i.e. all corrections applied).
%   domain - The geographic domain
%       It can be either a GeoDomain object, a string of the domain name, a cell array of the form {'domain name', buf},
%       or a N-by-2 matrix of longitude-latitude vertices.
%       The default value is 'oceans'.
%   L - The maximum degree of the Slepian expansion
%       The default value is 60.
%   timerange - The time range of the data
%       It can be either a 2-element datetime array or a 2-element numeric array.
%       The default value is [] (i.e. all time).
%
% Output arguments
%   slept - The Slepian expansions
%   dates - The dates of the data
%   G - The Slepian functions
%   V - The eigenvalues of the Slepian functions
%   N - The Shannon number
%   timeNum - The dates in datenum format
%
% Last modified by
%   2024/08/20, williameclee@arizona.edu (@williameclee)

function varargout = mascon2slept(varargin)
    %% Initialisation
    % Parse inputs
    [dataType, domain, L, timeRange, beQuiet] = parseinputs(varargin{:});

    %% Finding precomputed data
    outputFolder = fullfile(getenv('GRACEDATA'), 'SlepianExpansions');
    outputFile = sprintf('%s-CSRRL06%i-%s-SD-%s.mat', mfilename, L, domain.Id, dataType);
    outputPath = fullfile(outputFolder, outputFile);

    if exist(outputPath, 'file')

        load(outputPath, 'dataSlept', 'time');

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), outputPath);
        end

        if ~isempty(timeRange)
            timeRangeId = find(time >= timeRange(1) & time <= timeRange(2));
            time = time(timeRangeId);
            dataSlept = dataSlept(timeRangeId, :, :);
        end

        fakeData = zeros([addmup(L), 4]);
        [fakeData(:, 2), fakeData(:, 1)] = addmon(L);
        [dataSlept(1, :), V, N, ~, G] = plm2slep_new(fakeData, domain, L, "BeQuiet", true);
        timeNum = datenum(time); %#ok<DATNM>

        varargout = {dataSlept, time, G, V, N, timeNum};
        return;
    end

    %% Loading data
    switch dataType
        case 'deg123'
            [deg1Plmt, time] = mascon2plmt('deg1', L);
            masconC20Plmt = mascon2plmt('masconC20', L);
            masconC30Plmt = mascon2plmt('masconC30', L);
            slrC20Plmt = mascon2plmt('slrC20', L);
            [slrC30Plmt, timeSlrC30] = mascon2plmt('slrC30', L);

            dataPlmt = deg1Plmt;
            dataPlmt(:, :, 3:4) = dataPlmt(:, :, 3:4) ...
                - masconC20Plmt(:, :, 3:4) + slrC20Plmt(:, :, 3:4) ...
                - masconC30Plmt(:, :, 3:4);

            slrC30StartId = find(time == timeSlrC30(1), 1, 'first');
            dataPlmt(slrC30StartId:end, :, 3:4) = ...
                dataPlmt(slrC30StartId:end, :, 3:4) + slrC30Plmt(:, :, 3:4);
        otherwise
            [dataPlmt, time] = mascon2plmt(dataType, L);
    end

    timeNum = datenum(time); %#ok<DATNM>

    %% Computing Slepian expansions
    dataSlept = zeros([length(time), (L + 1) ^ 2]);

    [dataSlept(1, :), V, N, ~, G] = plm2slep_new(squeeze(dataPlmt(1, :, :)), domain, L, "BeQuiet", true);

    parfor t = 2:length(time)
        dataSlept(t, :) = plm2slep_new(squeeze(dataPlmt(t, :, :)), domain, L, "BeQuiet", true);
    end

    %% Saving and collecting output
    save(outputPath, 'dataSlept', 'time');

    if ~beQuiet
        fprintf('%s saved %s\n', upper(mfilename), outputPath);
    end

    if ~isempty(timeRange)
        timeRangeId = find(time >= timeRange(1) & time <= timeRange(2));
        time = time(timeRangeId);
        dataSlept = dataSlept(timeRangeId, :, :);
    end

    varargout = {dataSlept, time, G, V, N, timeNum};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addOptional(p, 'DataType', 'mascon', @ischar);
    addOptional(p, 'Domain', 'oceans', ...
        @(x) ischar(x) || isstring(x) || iscell(x) || ...
        isa(x, 'GeoDomain') || isnumeric(x));
    addOptional(p, 'L', 60, @isnumeric);
    addOptional(p, 'TimeRange', [], ...
        @(x) ((isdatetime(x) || isnumeric(x)) && length(x) == 2) || ...
        isempty(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));

    parse(p, varargin{:});
    dataType = p.Results.DataType;
    domain = p.Results.Domain;
    L = p.Results.L;
    timeRange = p.Results.TimeRange;
    beQuiet = p.Results.BeQuiet;

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    if isnumeric(timeRange)
        timeRange = datetime(timeRange, 'ConvertFrom', 'datenum');
    end

    varargout = {dataType, domain, L, timeRange, beQuiet};
end
