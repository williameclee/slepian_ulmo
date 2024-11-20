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
%   rslPlmt - Relative sea level spherical harmonic coefficients
%   time - Time stamps of the data, in datetime format
%
% See also
%   SOLVESLE
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)

function [rslPlmt, time] = grace2fingerprint(varargin)
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

    doGiaFeedback = ip.Results.GIAFeedback;
    doRotationFeedback = ip.Results.RotationFeedback;
    frame = ip.Results.Frame;

    beQuiet = ip.Results.BeQuiet;
    forceNew = ip.Results.ForceNew;
    saveData = ip.Results.SaveData;

    %% Data locating
    [filepath, fileexists] = findoutputpath(product, L, oceanDomain, frame, giaModel, doGiaFeedback, doRotationFeedback);

    if fileexists && ~forceNew
        load(filepath, 'rslPlmt', 'time');

        if ~isempty(Loutput) && Loutput < L
            rslPlmt = rslPlmt(1:addmup(Loutput), :, :);
        end

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), filepath);
        end

        return
    end

    %% Loading data
    [gracePlmt, ~, time] = grace2plmt_new(product{:}, 'SD', "BeQuiet", beQuiet);

    % Permute to [lmcosi, time]
    if size(gracePlmt, 1) == length(time)
        gracePlmt = permute(gracePlmt, [2, 3, 1]);
    end

    % Ignore uncertainty
    gracePlmt = gracePlmt(:, 1:4, :);

    % Convert time to datetime
    if isnumeric(time)
        time = datetime(time, "ConvertFrom", 'datenum');
    end

    % Load GIA model
    if doGia
        giaPlmt = gia2plmt(time, giaModel, L, "BeQuiet", beQuiet);

        % Permute to [lmcosi, time]
        if size(giaPlmt, 1) == length(time)
            giaPlmt = permute(giaPlmt, [2, 3, 1]);
        end

    end

    %% Data preprocessing
    % Make to the same degree
    if product{3} < L
        [order, degree] = addmon(L);
        gracePlmt(1:addmup(L), 1, :) = repmat(degree, [1, 1, length(time)]);
        gracePlmt(1:addmup(L), 2, :) = repmat(order, [1, 1, length(time)]);
    end

    % Correct for GIA
    gracePlmt(:, 3:4, :) = gracePlmt(:, 3:4, :) - giaPlmt(:, 3:4, :);
    [landloadingPlmt, Kocean] = localise(gracePlmt, oceanDomain, L, "Inverse", true);
    Kocean = eye(size(Kocean)) - Kocean;

    if doGiaFeedback
        landloadingPlmt(:, 3:4, :) = landloadingPlmt(:, 3:4, :) + giaPlmt(:, 3:4, :);
    end

    rslPlmt = nan(size(gracePlmt));

    parfor iTime = 1:length(time)
        rslPlmt(:, :, iTime) = solvesle(squeeze(landloadingPlmt(:, :, iTime)), L, oceanDomain, "Frame", frame, "RotationFeedback", doRotationFeedback, ...
            "OceanKernel", Kocean, "BeQuiet", true);
    end

    %% Saving and returning data
    if saveData
        save(filepath, 'rslPlmt', 'time');

        if ~beQuiet
            fprintf('%s saved %s\n', upper(mfilename), filepath);
        end

    end

    if ~isempty(Loutput) && Loutput < L
        rslPlmt = rslPlmt(1:addmup(Loutput), :, :);
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
