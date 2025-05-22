%% GRACE2SLEPT
% Takes GRACE/GRACE-FO gravimetry data created by GRACE2PLM and projects
% this data into the requested Slepian basis.
%
% Syntax
%   [date, slept] = grace2slept(product, domain, buf, L)
%   [date, slept] = ...
%       grace2slept(product, r, buf, L, phi, theta, omega, J)
%   [date, slept] = grace2slept(__, 'Name', Value)
%   [__, domain, G, eigfun, V, N] = grace2slept(__)
%
% Input arguments
%   product - The name of the data product
%       This is a cell array with three parts:
%       - the data center ('CSR', 'GFZ', or 'JPL'),
%       - the release level ('RL04', 'RL05', or 'RL06'),
%       - the dataproduct bandwidth.
%       The default product is {'CSR', 'RL06', 60}.
%   r - Radius of the concentration region in degrees
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A cell array of the domain of interest (string) and the buffer
%           value (scalar): {domain, buf}.
%       - An N-by-2 matrix of the domain of interest in [lon, lat] format.
%       - A GeoDomain object.
%   buf - Distance in degrees that the region outline will be enlarged
%              by BUFFERM [default: 0]
%   L - Bandwidth of the window [default: bandwidth of the data],
%              or bandpass (two degrees)
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the center of the tapers in degrees
%   J - Number of largest eigenfunctions in which to expand
%       - 'N': The Shannon number of eigenfunctions.
%       The default is all of them.
%   unit - The unit of the output data
%       - 'POT': Geopotential.
%       - 'SD': Surface mass density.
%       The default unit is 'POT'.
%   Deg1Correction, C20Correction, C30Correction - Wether you want to apply
%       the degree 1, C20, or C30 corrections
%       The default options are true.
%   TimeRange - The time range of GRACE data
%       A 2-by-1 datetime or timenum array specifying the starting and ending time of the GRACE data used.
%       The default range uses all data available.
%   ForceNew - Wether or not you want to force recomputation
%       The default option is false.
%   SaveData - Wether or not you want to save the output data
%       The default option is true.
%   BeQuiet - Wether or not you want to suppress output
%       The default option is false.
%
% Output and saved arguments
%   slept - The expansion coefficients of the GRACE data in specified unit
%       into the Slepian basispotential Slepian coefficients
%       The size is [nmonths x addmoff(Ldata)].
%   slept_err - The expansion coefficients of the calibrated errors into
%       the Slepian basis calibrated errors.
%       The size is [nmonths x addmoff(Ldata)].
%       For release RL06, the calibrated errors are not available.
%   date - Time stamps of the time series in DATETIME format
%   domain - The domain
%       If there was buffering, this will be a XY array of coordinates,
%       which you can use with SPHAREA to get the Shannon number.
%   G - The unitary matrix of localisation coefficients
%   eigfun - A cell array with cosine/sine coefficients eigenfunctions
%   V - The eigenvalues in this ordering
%   N - The Shannon number
%
% See also
%   GRACE2PLMT (GRACE2PLMT_NEW), PLM2SLEP
%
% Last modified by
%   2025/05/22, williameclee@arizona.edu (@williameclee)
%   2024/08/30, williameclee@arizona.edu (@williameclee)
%   2022/05/18, charig@princeton.edu (@harig00)
%   2012/06/26, fjsimons@alum.mit.edu (@fjsimons)

function varargout = grace2slept_new(varargin)
    %% Initialisation
    % Parse inputs
    [~, domain, L, phi, theta, omega, unit, timeRange, ...
         dataCentre, releaseLevel, Ldata, productStr, truncation, ...
         deg1corr, c20corr, c30corr, forceNew, saveData, beQuiet] = ...
        parseinputs(varargin{:});

    % Figure out if it's low-pass or band-pass
    bp = length(L) == 2;
    maxL = max(L);
    % The spherical harmonic dimension
    if ~bp
        ldim = (L + 1) ^ 2;
    else
        ldim = (L(2) + 1) ^ 2 - L(1) ^ 2;
    end

    % Check if you want the Shannon number of eigenfunctions
    if strcmp(truncation, 'N')
        truncation = ceil((L(end) + 1) ^ 2 * domain.Spharea);
    else
        truncation = conddefval(truncation, ldim);
    end

    % Output file
    [outputPath, outputExists] = getoutputfile(domain, L, ...
        productStr, truncation, unit, bp, ...
        deg1corr, c20corr, c30corr);

    %% Constructing the Slepian basis
    [~, ~, ~, lmcosiW, ~, ~, ~, ~, ~, ronmW] = addmon(maxL);

    if phi == 0 && theta == 0 && omega == 0
        [G, V, ~, ~, N] = glmalpha_new( ...
            domain, L, "J", truncation, "BeQuiet", beQuiet);
    else
        % Need to get a complete GLMALPHA but for the rotated basis
        % Definitely, "single-order" has lost its meaning here, but the MTAP
        % will still identify what the order of the unrotated original was
        [G, V, ~, ~, N] ...
            = glmalphapto(domain, L, phi, theta, omega);
        % Since GLMALPHAPTO currently has no option to limit a basis to J,
        % do it here
        G = G(:, 1:truncation);
    end

    % Sort by decreasing eigenvalue
    [V, sortId] = sort(V, 'descend');
    G = G(:, sortId);

    % If you don't do this, the eigenfunctions are ordered in the way
    %   that they correspond to single-orders back when, unrotated, they
    %   belonged to a polar cap, and the eigenvalues are sorted within
    %   these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.
    % Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
    CC = cell(1, size(G, 2));

    for j = 1:size(G, 2)
        % Create the blanks
        cosi = lmcosiW(:, 3:4);
        % Stick in the coefficients of the 1st eigentaper
        cosi(ronmW) = G(:, j);
        % Construct the full matrix
        CC{j} = [lmcosiW(:, 1:2) cosi];
    end

    % INITILIZATION COMPLETE

    % If this expansion already exists, load it.  Otherwise, or if we force
    % it, make a new one (e.g. if you added extra months to the database).
    if outputExists && ~forceNew
        load(outputPath, 'slept', 'dates', 'thedates', 'slepcoffs')

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), outputPath)
        end

        if ~exist('slept', 'var') && exist('slepcoffs', 'var')
            slept = slepcoffs;
        end

        if ~exist('dates', 'var') && exist('thedates', 'var')
            dates = thedates;
        end

        save(outputPath, 'slept', 'dates');

        % Truncate to required time range
        if ~isempty(timeRange)
            timeStartId = find(dates >= timeRange(1), 1, 'first');
            timeEndId = find(dates <= timeRange(2), 1, 'last');
            dates = dates(timeStartId:timeEndId);
            slept = slept(timeStartId:timeEndId, :);
        end

        dates = dates(:);
        dates = datetime(dates, 'ConvertFrom', 'datenum');
        varargout = {dates, slept, domain, G, CC, V, N};

        return
    end

    % Use GRACE2PLMT to get the GRACE data
    [sphCoeffs, ~, dates] = grace2plmt_new(dataCentre, releaseLevel, ...
        Ldata, unit, 'Deg1Correction', deg1corr, ...
        'C20Correction', c20corr, 'C30Correction', c30corr, "BeQuiet", beQuiet);
    dates = dates(:);
    % *** Here I changed this. Run grace2plmt once to update your data, and
    % then when you call forcenew=1 from now on it will just update the
    % expansion

    % Initialize new coefficients
    nDates = length(dates);
    slept = nan([nDates, truncation]);
    % slepcalerrors = nan(nmonths, truncation);

    % Limit everything to the window bandwidth
    potcoffsW = sphCoeffs(:, 1:addmup(L), 1:4);
    %cal_errorsW = cal_errors(:,1:size(lmcosiW,1),1:4);

    % Loop over the months
    for iDate = 1:nDates
        sphCoeffs = squeeze(potcoffsW(iDate, :, :));
        slept(iDate, :) = ...
            sphCoeffs(2 * size(sphCoeffs, 1) + ...
            ronmW(1:(maxL + 1) ^ 2))' * G;

        % Expand this month of CALIBRATED ERRORS into the Slepian basis
        %calerrors_month=squeeze(cal_errorsW(index,:,:));
        %slepcalerrors(index,:) = ...
        %    calerrors_month(2*size(calerrors_month,1)+ronmW(1:(maxL+1)^2))'*G;
    end

    if saveData
        save(outputPath, 'slept', 'dates');

        if ~beQuiet
            fprintf('%s saved %s\n', upper(mfilename), outputPath)
        end

    end

    % Truncate to required time range
    if ~isempty(timeRange)
        timeStartId = find(dates >= timeRange(1), 1, 'first');
        timeEndId = find(dates <= timeRange(2), 1, 'last');
        dates = dates(timeStartId:timeEndId);
        slept = slept(timeStartId:timeEndId, :);
    end

    % Collect output
    dates = dates(:);
    dates = datetime(dates, 'ConvertFrom', 'datenum');
    varargout = {dates, slept, domain, G, CC, V, N};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    productD = {'CSR', 'RL06', 60};
    domainD = 'greenland';
    LD = 18;
    taperD = 0;
    JD = [];
    unitD = 'SD';
    forceNewD = false;

    p = inputParser;
    addOptional(p, 'DataProduct', productD, ...
        @(x) (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'Domain', domainD, ...
        @(x) (ischar(x)) || isstring(x) || iscell(x) || ...
        isa(x, 'GeoDomain') || (isnumeric(x) && size(x, 2) == 2) || ...
        (isempty(x)));
    addOptional(p, 'L', LD, ...
        @(x) (isnumeric(x) && (length(x) <= 2)) || (isempty(x)));
    addOptional(p, 'phi', taperD, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'theta', taperD, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'omega', taperD, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'Truncation', JD, ...
        @(x) ((isnumeric(x) && isscalar(x) && x > 0) || ...
        strcmp(x, 'N')) || (isempty(x)));
    addOptional(p, 'Unit', unitD, ...
        @(x) (ischar(validatestring(x, {'POT', 'SD'}))) || (isempty(x)));
    addOptional(p, 'TimeRange', [], ...
        @(x) isempty(x) || ((isdatetime(x) || isnumeric(x))));
    addOptional(p, 'ForceNew', forceNewD, ...
        @(x) (isnumeric(x) && (x == 0 || x == 1)) || islogical(x) ...
        || (isempty(x)));
    addOptional(p, 'MoreRegionSpecs', {}, @iscell);
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveData', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'Deg1Correction', true, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'C20Correction', true, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'C30Correction', true, ...
        @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    product = conddefval(p.Results.DataProduct, productD);
    domain = conddefval(p.Results.Domain, domainD);

    L = conddefval(p.Results.L, LD);
    phi = conddefval(p.Results.phi, taperD);
    theta = conddefval(p.Results.theta, taperD);
    omega = conddefval(p.Results.omega, taperD);
    unit = conddefval(p.Results.Unit, unitD);
    domainSpecs = p.Results.MoreRegionSpecs;
    J = conddefval(p.Results.Truncation, JD);
    forceNew = conddefval(logical(p.Results.ForceNew), forceNewD);
    beQuiet = logical(p.Results.BeQuiet);
    saveData = logical(p.Results.SaveData);
    deg1correction = logical(p.Results.Deg1Correction);
    c20correction = logical(p.Results.C20Correction);
    c30correction = logical(p.Results.C30Correction);
    timeRange = p.Results.TimeRange;

    if isnumeric(product{2})
        product{2} = ['RL0', num2str(product{2})];
    end

    dataCentre = product{1};
    releaseLevel = product{2};
    Ldata = product{3};
    productId = [dataCentre, releaseLevel, num2str(Ldata)];

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, domainSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2}, ...
            domainSpecs{:});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:}, domainSpecs{:});
    end

    if isdatetime(timeRange)
        timeRange = datenum(timeRange); %#ok<DATNM>
    end

    varargout = ...
        {product, domain, L, ...
         phi, theta, omega, unit, timeRange, ...
         dataCentre, releaseLevel, Ldata, productId, J, ...
         deg1correction, c20correction, c30correction, ...
         forceNew, saveData, beQuiet};
end

function [outputPath, outputExists] = getoutputfile(domain, L, ...
        productId, truncation, unit, bp, ...
        deg1corr, c20corr, c30corr)

    % Folder
    if ~isempty(getenv('GRACE'))
        outputFolder = fullfile(getenv('GRACE'), ...
        'SlepianExpansions');
    else
        outputFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'SlepianExpansions');
    end

    % File name
    if bp
        Lstr = sprintf('%i-%i', L(1), L(2));
    else
        Lstr = sprintf('%i', L);
    end

    switch domainType(domain)
        case 'polar'
            r = domain;

            outputFile = sprintf( ...
                'grace2slept-%s-CAP-%i-%s-%i-%s.mat', ...
                productId, r, Lstr, truncation, unit);

        case {'geodomain', 'lonlat'}

            switch domainType(domain)
                case 'geodomain'
                    domainId = domain.Id;
                case 'lonlat'
                    domainId = hash(domain, 'sha1');
            end

            % The name of the save file
            outputFile = sprintf( ...
                'grace2slept-%s-%s-%s-%i-%s.mat', ...
                productId, domainId, Lstr, truncation, unit);

    end

    if ~deg1corr
        outputFile = strrep(outputFile, '.mat', '-nDeg1.mat');
    end

    if ~c20corr
        outputFile = strrep(outputFile, '.mat', '-nC20.mat');
    end

    if ~c30corr
        outputFile = strrep(outputFile, '.mat', '-nC30.mat');
    end

    outputPath = fullfile(outputFolder, outputFile);

    outputExists = isfile(outputPath);

end

function domainType = domainType(domain)

    if isa(domain, 'GeoDomain')
        domainType = 'geodomain';
    elseif isnumeric(domain) && isscalar(domain)
        domainType = 'polar';
    elseif isnumeric(domain) && size(domain, 2) == 2
        domainType = 'lonlat';
    else
        error('Unknown domain type')
    end

end
