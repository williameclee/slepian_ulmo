%% SLF2SLEPT
% Computes Slepian expansions of the sea-level fingerprints (SLF)
%
% Syntax
%   [time, sleptLandload, sleptRsl, sleptGeoid, sleptBedrock] = ...
% 		slf2slept(dataType, domain, L)
%   [__] = slf2slept(__, timerange)
%   [__] = slf2slept(__, 'Name', Value)
%
% Input arguments
%   dataType - The type of data to use
%       It consists of three elements:
%       - Data centre: 'CSR', 'JPL', or 'GFZ'
%       - Reference frame: 'CF' (centre-of-figure) or 'CM' (centre-of-mass)
%       - Rotation: true or false
%		This can be speficied as three arguments or a cell array.
%       The default data type is {'CSR', 'CM', true}.
%   domain - The geographic domain
%       It can be either a GeoDomain object, a string of the domain name, a cell array of the form {'domain name', buf},
%       or a N-by-2 matrix of longitude-latitude vertices.
%       The default value is 'oceans'.
%   L - The maximum degree of the Slepian expansion
%       The default value is 60.
%   timerange - The time range of the data
%       It can be either a 2-element datetime array or a 2-element numeric array.
%       The default value is [] (i.e. all time).
%	ForceNew - Whether to force the recomputation of the Slepian expansions
%		The default option is false.
%	SaveData - Whether to save the Slepian expansions
%		The default option is true.
%	BeQuiet - Whether to suppress the output
%		The default option is soft quiet (i.e. only surpress outputs from other functions).
%
% Output arguments
%   time - The dates of the data
%   sleptLandload, sleptRsl, sleptGeoid, sleptBedrock - The Slepian expansions of different components of the SLF
%
% Last modified by
%   2024/08/22, williameclee@arizona.edu (@williameclee)

function varargout = slf2slept(varargin)
    %% Initialisation
    % Parse inputs
    [pcenter, rframe, rotation, domain, L, timeRange, beQuiet, forceNew, saveData] = parseinputs(varargin{:});

    %% Finding precomputed data
    outputFolder = fullfile(getenv('GRACEDATA'), 'SlepianExpansions');
    outputFile = sprintf('%s-SLF%i-%s.mat', upper(mfilename), L, domain.Id);
    outputPath = fullfile(outputFolder, outputFile);

    if exist(outputPath, 'file') && ~forceNew

        load(outputPath, 'time', 'sleptLandload', 'sleptRsl', 'sleptGeoid', 'sleptSd', 'sleptBedrock');

        if beQuiet <= 2
            fprintf('%s loaded %s\n', upper(mfilename), outputPath);
        end

    else

        %% Loading data
        [time, plmtLandload, plmtRsl, plmtGeoid, plmtSd, plmtBedrock] = ...
            slf2plmt(pcenter, rframe, rotation, L);

        %% Computing Slepian expansions
        % Preallocate
        sleptLandload = zeros([length(time), (L + 1) ^ 2]);
        sleptRsl = zeros([length(time), (L + 1) ^ 2]);
        sleptGeoid = zeros([length(time), (L + 1) ^ 2]);
        sleptSd = zeros([length(time), (L + 1) ^ 2]);
        sleptBedrock = zeros([length(time), (L + 1) ^ 2]);

        % Compute Slepian expansions
        parfor t = 1:length(time)
            sleptLandload(t, :) = plm2slep_new(squeeze(plmtLandload(:, :, t)), domain, L, "BeQuiet", beQuiet);
            sleptRsl(t, :) = plm2slep_new(squeeze(plmtRsl(:, :, t)), domain, L, "BeQuiet", beQuiet);
            sleptGeoid(t, :) = plm2slep_new(squeeze(plmtGeoid(:, :, t)), domain, L, "BeQuiet", beQuiet);
            sleptSd(t, :) = plm2slep_new(squeeze(plmtSd(:, :, t)), domain, L, "BeQuiet", beQuiet);
            sleptBedrock(t, :) = plm2slep_new(squeeze(plmtBedrock(:, :, t)), domain, L, "BeQuiet", beQuiet);
        end

        if saveData
            save(outputPath, 'time', 'sleptLandload', 'sleptRsl', 'sleptGeoid', 'sleptSd', 'sleptBedrock', '-v7');

            if ~beQuiet
                fprintf('%s saved %s\n', upper(mfilename), outputPath);
            end

        end

    end

    %% Post-processing
    [~, sleptLandload, sleptRsl, sleptGeoid, sleptSd, sleptBedrock] = ...
        truncatetimerange(time, timeRange, sleptLandload, sleptRsl, sleptGeoid, sleptSd, sleptBedrock);

    %% Collecting and displaying outputs
    varargout = {time, sleptLandload, sleptRsl, sleptGeoid, sleptSd, sleptBedrock};

end

%% Subfunctions
function varargout = parseinputs(varargin)

    if isempty(varargin)
        pcenter = 'CSR';
        rframe = 'CM';
        rotation = true;
    else
        productParser = inputParser;
        addOptional(productParser, 'ProcudtCenter', 'CSR', ...
            @(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'CSR', 'JPL', 'GFZ'}));
        addOptional(productParser, 'ReferenceFrame', 'CM', ...
            @(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'CF', 'CM'}));
        addOptional(productParser, 'Rotation', true, ...
            @(x) islogical(x) || isnumeric(x));

        if iscell(varargin) && length(varargin{1}) == 3
            parse(productParser, varargin{1}{:});
            varargin = varargin(2:end);
        else
            parse(productParser, varargin{1:3});
            varargin = varargin(4:end);
        end

        pcenter = upper(productParser.Results.ProcudtCenter);
        rframe = upper(productParser.Results.ReferenceFrame);
        rotation = logical(productParser.Results.Rotation);
    end

    if ~ismember(pcenter, {'CSR', 'JPL', 'GFZ'})
        error('Unknown pcenter: %s, must be CSR, JPL, or GFZ', pcenter);
    end

    if any(strcmpi(rframe, {'FIGURE', 'COF'}))
        rframe = 'CF';
    elseif any(strcmpi(rframe, {'MASS', 'COM'}))
        rframe = 'CM';
    elseif ~any(strcmpi(rframe, {'CF', 'CM'}))
        error('Unknown reference frame: %s, must be CF or CM', rframe);
    end

    p = inputParser;
    addOptional(p, 'Domain', 'oceans', ...
        @(x) ischar(x) || isstring(x) || iscell(x) || ...
        isa(x, 'GeoDomain') || isnumeric(x));
    addOptional(p, 'L', 60, @isnumeric);
    addOptional(p, 'TimeRange', [], ...
        @(x) ((isdatetime(x) || isnumeric(x)) && length(x) == 2) || ...
        isempty(x));
    addParameter(p, 'BeQuiet', 0.5, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveData', true, @(x) islogical(x) || isnumeric(x));

    parse(p, varargin{:});
    domain = p.Results.Domain;
    L = p.Results.L;
    timeRange = p.Results.TimeRange;
    beQuiet = uint8(2 * p.Results.BeQuiet);
    forceNew = p.Results.ForceNew;
    saveData = p.Results.SaveData;

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

    forceNew = uint8(forceNew);

    varargout = {pcenter, rframe, rotation, domain, L, timeRange, beQuiet, forceNew, saveData};
end

function varargout = truncatetimerange(time, trange, varargin)

    if isempty(time)
        varargout = [{}, varargin];
        return
    end

    maxTrange = [min(time), max(time)];

    if isscalar(trange)
    elseif isempty(trange) || all(isnat(trange))
        trange = maxTrange;
    elseif any(isnat(trange)) || ...
            trange(1) < maxTrange(1) || trange(2) > maxTrange(2)
        trange = ...
            [max(trange(1), maxTrange(1), "omitmissing"), ...
             min(trange(2), maxTrange(2), "omitmissing")];
        warning('The timerange is partially invalid or outside the range of the data files, truncated to match the data files:\n%s - %s', ...
            datetime(trange(1), "Format", 'yyyy/MM/dd'), ...
            datetime(trange(2), "Format", 'yyyy/MM/dd'));
    end

    if isscalar(trange)
        %#ok<*DATNM>
        timeId = find(min(datenum(time(:) - trange), [], 'omitnan'), 1);
    else
        timeId = time >= trange(1) & time <= trange(2);
    end

    time = time(timeId);

    if nargin == 2
        varargout = {time};
        return
    end

    varargin = varargin(:)';

    for i = 1:length(varargin)
        varargin{i} = squeeze(varargin{i}(timeId, :));
    end

    varargout = [{time}, varargin];
end
