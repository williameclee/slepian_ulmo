%% GRACE2TREND
% Reads GRACE data and computes the mass trend projected onto a set of Slepian functions.
%
% Last modified by
%   2024/08/30, williameclee@arizona.edu (@williameclee)

function varargout = grace2trend(domain, varargin)
    %% Initialisation
    [domain, L, timeRange, fitwhat, product, gia, deg1cor, c20cor, ...
         c30cor, ~, beQuiet] = parseinputs(domain, varargin{:});

    %% Loading data
    [date, slept, ~, G, ~, V, N] = ...
        grace2slept_new(product, domain, L, ...
        "TimeRange", timeRange, "Unit", 'SD', ...
        "Deg1Correction", deg1cor, "C20Correction", c20cor, ...
        "C30Correction", c30cor, "BeQuiet", beQuiet);

    % Correct for GIA
    if gia
        [~, gisSlept] = ...
            gia2slept(date, gia, domain, L, "BeQuiet", beQuiet);
        slept = slept - gisSlept;
    end

    %% Computing trends
    [sleptSig, sleptRes, ~, ~, total, alphavarall, fitparams, ...
         fitparamerrors, totalfit] = ...
        slept2resid_new(slept, date, fitwhat, "Domain", domain, ...
        "Unit", 'year', "BeQuiet", beQuiet);

    totalerror = sqrt(alphavarall);

    varargout = ...
        {date, total, totalerror, totalfit, fitparams, fitparamerrors, ...
         slept, sleptSig, sleptRes, N, G, V};
end

%% Subfunctions
function varargout = parseinputs(varargin)
    fitwhatD = [3, 365.25, 182.625, 161];
    p = inputParser;
    addRequired(p, 'Domain', ...
        @(x) ischar(x) || isstring(x) || iscell(x) || isa(x, 'GeoDomain'));
    addOptional(p, 'L', 30, @isnumeric);
    addOptional(p, 'TimeRange', ...
        [datetime(2003, 1, 1), datetime(2022, 12, 31)], ...
        @(x) isdatetime(x) && length(x) == 2);
    addOptional(p, 'FitWhat', fitwhatD, ...
        @(x) isnumeric(x));
    addOptional(p, 'DataProduct', {'CSR', 'RL06', 60}, ...
        @(x) iscell(x) && length(x) == 3 && ...
        ischar(x{1}) && ischar(x{2}) && isnumeric(x{3}));
    addParameter(p, 'GIACorrection', 'Paulson07', ...
        @(x) ischar(x) || islogical(x));
    addParameter(p, 'Deg1Correction', true, @islogical);
    addParameter(p, 'C20Correction', true, @islogical);
    addParameter(p, 'C30Correction', true, @islogical);
    addParameter(p, 'Truncation', [], @isnumeric);
    addParameter(p, "BeQuiet", 0.5, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    domain = p.Results.Domain;
    L = p.Results.L;
    timeRange = p.Results.TimeRange;
    fitwhat = conddefval(p.Results.FitWhat, fitwhatD);
    DataProduct = p.Results.DataProduct;
    giaModel = p.Results.GIACorrection;
    deg1cor = p.Results.Deg1Correction;
    c20cor = p.Results.C20Correction;
    c30cor = p.Results.C30Correction;
    truncation = p.Results.Truncation;
    beQuiet = uint8((double(p.Results.BeQuiet) * 2));

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
        timeRange = datetime(timeRange, "ConvertFrom", 'datenum');
    end

    if islogical(giaModel) && giaModel
        giaModel = 'Paulson07';
    end

    varargout = ...
        {domain, L, timeRange, fitwhat, DataProduct, giaModel, ...
         deg1cor, c20cor, c30cor, truncation, beQuiet};
end
