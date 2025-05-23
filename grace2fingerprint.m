%% GRACE2FINGERPRINT
% Computes the relative sea level fingerprint based on forcing from GRACE data.
%
% Syntax
%   [rslPlmt, time] = grace2fingerprint(L)
%   [rslPlmt, time] = grace2fingerprint(L, "GIA", model, "GIAFeedback", feedback)
%   [rslPlmt, time] = grace2fingerprint(__, "RotationFeedback", feedback)
%   [rslPlmt, time] = grace2fingerprint(__, "Name", Value)
%
% Input arguments
%   product - GRACE product
%       A cell array of the form {centre, release, degree}.
%       The default product is {'CSR', 'RL06', 60}
%   L - Bandwidth to compute spherical harmonic coefficients
%       The default value is 96
%   oceanDomain - Ocean domain
%       A geographic domain (GeoDomain object).
%       The default domain is all oceans with a buffer of 1 degree
%   "GIA" - Glacial isostatic adjustment model
%       The default model is 'ICE-6G_D'
%   "GIAFeedback" - Whether to include GIA feedback
%       The default value is true
%   "RotationFeedback" - Whether to include rotation feedback
%       The default value is true
%   "Frame" - Reference frame
%       Centre of mass (CM) frame or centre of figure (CF) frame.
%       The default frame is the CF frame, which is the frame used in slepian_delta
%   "Truncation" - Truncate the output to a lower degree
%       The default value is the same as L
%
% Output arguments
%   rslPlmt - Relative sea level load spherical harmonic coefficients
%   time - Time stamps of the data, in datetime format
%
% Notes
%   It doesn't make much sense to compute the GIA feedback or steric correction at this moment, so both arguments are disabled.
%
% See also
%   SOLVESLE
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/03/28, williameclee@arizona.edu (@williameclee)

function [rslLoadSpht, time] = grace2fingerprint(varargin)
    %% Parsing inputs
    ip = inputParser;
    addOptional(ip, 'Product', {'CSR', 'RL06', 60}, ...
        @(x) iscell(x) && length(x) == 3 && all(cellfun(@ischar, x(1:2))) && isnumeric(x{3}));
    addOptional(ip, 'L', 96, @(x) isscalar(x) && isnumeric(x));
    addOptional(ip, 'OceanDomain', GeoDomain('alloceans', "Buffer", 0.5), ...
        @(x) isa(x, 'GeoDomain') || ischar(x) || (iscell(x) && length(x) == 2) || (isnumeric(x) && size(x, 2) == 2));
    addOptional(ip, 'GIA', 'ice6gd', @(x) ischar(x) || islogical(x));
    addParameter(ip, 'GIAFeedback', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
    addParameter(ip, 'RotationFeedback', true, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
    % addParameter(ip, 'StericCorrection', false, ...
    %     @(x) (iscell(x) && length(x) == 2 && all(cellfun(@ischar, x))) || islogical(x));
    addParameter(ip, 'RecomputeDegree1', false, @(x) islogical(x));
    addParameter(ip, 'Frame', 'CF', @(x) ischar(validatestring(x, {'CF', 'CM'})));
    addParameter(ip, 'Truncation', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'TimeRange', [], @(x) (isdatetime(x) || isnumeric(x)) && length(x) == 2 || isempty(x));
    addParameter(ip, 'BeQuiet', false, @(x) isscalar(x));
    addParameter(ip, 'ForceNew', false, @(x) isscalar(x));
    addParameter(ip, 'SaveData', true, @(x) isscalar(x));
    parse(ip, varargin{:});

    beQuiet = logical(ip.Results.BeQuiet);
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    product = ip.Results.Product;
    L = ip.Results.L;
    Loutput = ip.Results.Truncation;
    oceanDomain = ip.Results.OceanDomain;
    giaModel = lower(ip.Results.GIA);
    % stericModel = ip.Results.StericCorrection;
    redoDeg1 = logical(ip.Results.RecomputeDegree1);

    if islogical(giaModel) && ~giaModel
        doGia = false;
    else
        doGia = true;

        if islogical(giaModel)
            giaModel = 'ice6gd';

            if ~beQuiet
                fprintf('%s using default GIA model %s\n', upper(mfilename), giaModel);
            end

        end

    end

    % if islogical(stericModel) && ~stericModel
    %     doSteric = false;
    % else
    %     doSteric = true;

    %     if islogical(stericModel)
    %         stericModel = {'NOAA', '2000total'};

    %         if ~beQuiet
    %             fprintf('%s using default steric model %s\n', upper(mfilename), strjoin(stericModel, '_'));
    %         end

    %     end

    % end

    timelim = ip.Results.TimeRange;

    if isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    doGiaFeedback = logical(ip.Results.GIAFeedback);
    doRotationFeedback = logical(ip.Results.RotationFeedback);
    frame = ip.Results.Frame;

    %% Data locating
    [filepath, fileexists] = findoutputpath(product, L, oceanDomain, frame, giaModel, doGiaFeedback, doRotationFeedback, redoDeg1);
    % [filepath, fileexists] = findoutputpath(product, L, oceanDomain, frame, giaModel, doGiaFeedback, doRotationFeedback, stericModel);

    if fileexists && ~forceNew

        load(filepath, 'rslLoadTspht', 'time');

        if ~exist('rslLoadTspht', 'var')
            load(filepath, 'rslTrunPlmt', 'time');
            rslLoadTspht = rslTrunPlmt;
        end

        if ~isempty(timelim)
            rslLoadTspht = rslLoadTspht(:, :, time >= timelim(1) & time <= timelim(2));
            time = time(time >= timelim(1) & time <= timelim(2));
        end

        nData = size(rslLoadTspht, 3);
        [order, degree] = addmon(L);
        rslLoadSpht = cat(2, repmat(degree, [1, 1, nData]), repmat(order, [1, 1, nData]), rslLoadTspht);

        if ~isempty(Loutput) && Loutput < L
            rslLoadSpht = rslLoadSpht(1:addmup(Loutput), :, :);
        end

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), filepath);
        end

        return
    end

    %% Loading data
    deg1Args = false;

    if redoDeg1
        deg1Args = {60, L, giaModel, GeoDomain('alloceans', "Buffer", 1.5)};
    end

    [graceSpht, ~, time] = grace2plmt_new(product{:}, "RecomputeDegree1", deg1Args, ...
        "OutputFormat", 'traditional', "BeQuiet", beQuiet);

    % Make sure the format is [lmcosi, time]
    % if size(graceSpht, 1) == length(time)
    %     graceSpht = permute(graceSpht, [2, 3, 1]);
    % end

    % Ignore uncertainty
    graceSpht = graceSpht(:, 1:4, :);
    % Get the variation
    graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) - mean(graceSpht(:, 3:4, :), 3);

    % Convert time to datetime
    if isnumeric(time)
        time = datetime(time, "ConvertFrom", 'datenum');
    end

    % Load GIA model
    if doGia
        giaSdSpht = gia2plmt(time, giaModel, L, "OutputFormat", 'traditional', "BeQuiet", beQuiet);
    end

    % % Load steric model
    % if doSteric
    %     [stericSphtIn, stericTime, hasStericDataIn] = steric2plmt(stericModel{:}, "L", L, "TimeRange", [time(1), time(end)], 'BeQuiet', beQuiet);
    %     hasStericData = false(size(time));
    %     stericSpht = zeros([addmup(L), 4, length(time)]);

    %     for iTime = 1:length(time)
    %         iTimeIn = find(stericTime == time(iTime));

    %         if hasStericDataIn(iTimeIn)
    %             hasStericData(iTime) = true;
    %             stericSpht(:, :, iTime) = stericSphtIn(:, :, iTimeIn);
    %         end

    %     end

    % end

    %% Data preprocessing
    % Make to the same degree
    if product{3} < L
        [order, degree] = addmon(L);
        graceSpht(1:addmup(L), 1, :) = repmat(degree, [1, 1, length(time)]);
        graceSpht(1:addmup(L), 2, :) = repmat(order, [1, 1, length(time)]);
    end

    % Correct for GIA
    if doGia
        graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) - giaSdSpht(:, 3:4, :);
    end

    [forcingSpht, Kocean] = localise(graceSpht, oceanDomain, L, "Inverse", true);
    Kocean = eye(size(Kocean)) - Kocean;

    if doGiaFeedback
        giaGeoidSpht = gia2plmt(time, giaModel, L, "OutputFormat", 'traditional', "OutputField", 'geoid', "BeQuiet", beQuiet);
        giaVlmSpht = giaz2plmt(time, giaModel, L, "BeQuiet", beQuiet);
        giaVlmSpht(:, 3:4, :) = giaVlmSpht(:, 3:4, :) / 1000; % mm -> m

        % Get the variation
        giaGeoidSpht(:, 3:4, :) = giaGeoidSpht(:, 3:4, :) - mean(giaGeoidSpht(:, 3:4, :), 3);
        giaVlmSpht(:, 3:4, :) = giaVlmSpht(:, 3:4, :) - mean(giaVlmSpht(:, 3:4, :), 3);
    else
        giaGeoidSpht = double.empty([0, 0, length(time)]);
        giaVlmSpht = double.empty([0, 0, length(time)]);
    end

    % if doSteric
    %     stericSpht(:, 3:4, hasStericData) = stericSpht(:, 3:4, hasStericData) - mean(stericSpht(:, 3:4, hasStericData), 3);
    % else
    %     stericSpht = double.empty([0, 0, length(time)]);
    % end

    rslLoadSpht = zeros(size(graceSpht));

    parfor iTime = 1:size(rslLoadSpht, 3)
        rslLoadSpht(:, :, iTime) = ...
            solvesle(forcingSpht(:, :, iTime), L, oceanDomain, ...
            "Frame", frame, "RotationFeedback", doRotationFeedback, ...
            "OceanKernel", Kocean, "BeQuiet", true);
    end

    %% Saving and returning data
    if saveData
        rslLoadTspht = rslLoadSpht(:, 3:4, :);
        save(filepath, 'rslLoadTspht', 'time');

        if ~beQuiet
            fprintf('%s saved %s\n', upper(mfilename), filepath);
        end

    end

    if ~isempty(timelim)
        rslLoadSpht = rslLoadSpht(:, :, time >= timelim(1) & time <= timelim(2));
        time = time(time >= timelim(1) & time <= timelim(2));
    end

    if ~isempty(Loutput) && Loutput < L
        rslLoadSpht = rslLoadSpht(1:addmup(Loutput), :, :);
    end

end

%% Subfunctions
function [filepath, fileexists] = findoutputpath(product, L, oceanDomain, frame, giaModel, doGiaFeedback, doRotationFeedback, redoDeg1)
    productName = sprintf("%s_%s_%d", product{:});

    if isequal(oceanDomain, GeoDomain('alloceans', "Buffer", 0.5))
        oceanName = '';
    else

        switch class(oceanDomain)
            case 'GeoDomain'
                oceanName = sprintf('-%s', oceanDomain.Id);
            case 'char'
                oceanName = sprintf('-%s', capitalise(oceanDomain));
            case 'cell'
                oceanName = sprintf('-%s%d', capitalise(oceanDomain{1}), oceanDomain{2});
            case 'numeric'
                oceanName = sprintf('-%s', hash(oceanDomain, 'sha1'));
        end

    end

    if islogical(giaModel) && ~giaModel
        giaName = '';
    else
        giaName = sprintf('-%s', giaModel);
    end

    % if islogical(stericModel) && ~stericModel
    %     stericName = '';
    % else
    %     stericName = sprintf('-%s', strjoin(stericModel, '_'));
    % end

    if redoDeg1
        deg1Name = '-Rdeg1';
    else
        deg1Name = '';
    end

    filefolder = fullfile(getenv('IFILES'), 'FINGERPRINT');

    if ~exist(filefolder, 'dir')
        mkdir(filefolder);
    end

    filename = sprintf('%s%s-L%d%s-%s-G%d-R%d%s.mat', productName, giaName, L, oceanName, frame, doGiaFeedback, doRotationFeedback, deg1Name);
    filepath = fullfile(filefolder, filename);
    fileexists = isfile(filepath);
end
