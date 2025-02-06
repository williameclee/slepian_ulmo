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
% See also
%   SOLVESLE
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/02/05, williameclee@arizona.edu (@williameclee)

function [rslLoadSpht, time] = grace2fingerprint(varargin)
    %% Parsing inputs
    ip = inputParser;
    addOptional(ip, 'Product', {'CSR', 'RL06', 60}, ...
        @(x) iscell(x) && length(x) == 3 && all(cellfun(@ischar, x(1:2))) && isnumeric(x{3}));
    addOptional(ip, 'L', 96, @(x) isscalar(x) && isnumeric(x));
    addOptional(ip, 'OceanDomain', GeoDomain('alloceans', "Buffer", 1), ...
        @(x) isa(x, 'GeoDomain') || ischar(x) || (iscell(x) && length(x) == 2) || (isnumeric(x) && size(x, 2) == 2));
    addOptional(ip, 'GIA', 'ice6gd', @(x) ischar(x) || islogical(x));
    addParameter(ip, 'GIAFeedback', true, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
    addParameter(ip, 'RotationFeedback', true, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
    addParameter(ip, 'Frame', 'CF', @(x) ischar(validatestring(x, {'CF', 'CM'})));
    addParameter(ip, 'Truncation', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'TimeRange', [], @(x) (isdatetime(x) || isnumeric(x)) && length(x) == 2);
    addParameter(ip, 'BeQuiet', false, @(x) isscalar(x) && islogical(x));
    addParameter(ip, 'ForceNew', false, @(x) isscalar(x) && islogical(x));
    addParameter(ip, 'SaveData', true, @(x) isscalar(x) && islogical(x));
    parse(ip, varargin{:});
    product = ip.Results.Product;
    L = ip.Results.L;
    Loutput = ip.Results.Truncation;
    oceanDomain = ip.Results.OceanDomain;
    giaModel = lower(ip.Results.GIA);

    if islogical(giaModel) && ~giaModel
        doGia = false;
    else
        doGia = true;
    end

    timelim = ip.Results.TimeRange;

    if isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    doGiaFeedback = ip.Results.GIAFeedback;
    doRotationFeedback = ip.Results.RotationFeedback;
    frame = ip.Results.Frame;

    beQuiet = ip.Results.BeQuiet;
    forceNew = ip.Results.ForceNew;
    saveData = ip.Results.SaveData;

    %% Data locating
    [filepath, fileexists] = findoutputpath(product, L, oceanDomain, frame, giaModel, doGiaFeedback, doRotationFeedback);

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
    [graceSpht, ~, time] = grace2plmt_new(product{:}, 'SD', "BeQuiet", beQuiet);

    % Permute to [lmcosi, time]
    if size(graceSpht, 1) == length(time)
        graceSpht = permute(graceSpht, [2, 3, 1]);
    end

    % Ignore uncertainty
    graceSpht = graceSpht(:, 1:4, :);
    graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) - mean(graceSpht(:, 3:4, :), 3);

    % Convert time to datetime
    if isnumeric(time)
        time = datetime(time, "ConvertFrom", 'datenum');
    end

    % Load GIA model
    if doGia
        giaSdSpht = gia2plmt(time, giaModel, L, "OutputFormat", 'traditional', "BeQuiet", beQuiet);
    end

    %% Data preprocessing
    % Make to the same degree
    if product{3} < L
        [order, degree] = addmon(L);
        graceSpht(1:addmup(L), 1, :) = repmat(degree, [1, 1, length(time)]);
        graceSpht(1:addmup(L), 2, :) = repmat(order, [1, 1, length(time)]);
    end

    % Correct for GIA
    graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) - giaSdSpht(:, 3:4, :);
    [forcingSpht, Kocean] = localise(graceSpht, oceanDomain, L, "Inverse", true);
    Kocean = eye(size(Kocean)) - Kocean;

    if doGiaFeedback
        giaGeoidSpht = gia2plmt(time, giaModel, L, "OutputFormat", 'traditional', "OutputField", 'geoid', "BeQuiet", beQuiet);
        giaVlmSpht = giaz2plmt(time, giaModel, L, "BeQuiet", beQuiet);
        giaVlmSpht(:, 3:4, :) = giaVlmSpht(:, 3:4, :) / 1000; % mm -> m
    else
        giaGeoidSpht = [];
        giaVlmSpht = [];
    end

    rslLoadSpht = zeros(size(graceSpht));

    if doGiaFeedback

        parfor iTime = 1:size(rslLoadSpht, 3)
            rslLoadSpht(:, :, iTime) = ...
                solvesle(squeeze(forcingSpht(:, :, iTime)), L, oceanDomain, ...
                "Frame", frame, "RotationFeedback", doRotationFeedback, ...
                "GiaGeoidSph", squeeze(giaGeoidSpht(:, :, iTime)), "GiaVlmSph", squeeze(giaVlmSpht(:, :, iTime)), ...
                "OceanKernel", Kocean, "BeQuiet", true);

        end

    else

        parfor iTime = 1:size(rslLoadSpht, 3)
            rslLoadSpht(:, :, iTime) = ...
                solvesle(squeeze(forcingSpht(:, :, iTime)), L, oceanDomain, ...
                "Frame", frame, "RotationFeedback", doRotationFeedback, ...
                "OceanKernel", Kocean, "BeQuiet", true);

        end

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
function [filepath, fileexists] = findoutputpath(product, L, oceanDomain, frame, giaModel, doGiaFeedback, doRotationFeedback)
    productName = sprintf("%s_%s_%d", product{:});

    if isequal(oceanDomain, GeoDomain('alloceans', "Buffer", 1))
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

    filefolder = fullfile(getenv('IFILES'), 'FINGERPRINT');

    if ~exist(filefolder, 'dir')
        mkdir(filefolder);
    end

    filename = sprintf('%s%s-L%d%s-%s-G%d-R%d.mat', productName, giaName, L, oceanName, frame, doGiaFeedback, doRotationFeedback);
    filepath = fullfile(filefolder, filename);
    fileexists = isfile(filepath);
end
