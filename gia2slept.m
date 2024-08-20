%% GIA2SLEPT
% Computes the GIA correction time series projected onto the Slepian basis.
%
% Syntax
%   [time, slept] = gia2slept(model, domain)
%   [time, slept] = gia2slept(model, domain, L)
%   [time, slept] = gia2slept(model, r, L, phi, theta, omega)
%   [time, slept] = gia2slept(days, __)
%   [time, slept] = gia2slept(time, __)
%   [time, slept] = gia2slept(__, 'Name', Value)
%   [time, slept, sleptU, sleptL, total] = gia2slept(__)
%
% Input arguments
%   model - Name of the GIA model
%       - 'Steffen_ice6g_vm5a': A model computed by H. Steffen using the
%           ice6g ice model and vm5a viscosity profile. Other models from
%           this dataset are also available and use the original naming
%           scheme. For example, 'Steffen_anu-ice_i72'.
%           This family of models can also be specified as a cell array,
%           e.g. {'Steffen', 'ice6g', 'vm5a'}.
%       - 'Paulson07': A model based on the ICE-5G ice load model of
%           Peltier (2004). Suitable for both Antarctica and Greenland. As
%           corrected by Geruo A and J. Wahr.
%           Please avoid using this model for oceans.
%       - 'Wangetal08': A model based on the older ICE-4G ice model, and
%           viscosity which varies laterally. Suitable for Greenland.
%       - 'IJ05_R2': A model based on the Ivins et al. (2013) ice model.
%           Suitable for Antarctica.
%       - 'IJ05': A model based on the Ivins and James (2005) ice model.
%           Suitable for Antarctica.
%       - 'W12a_v1': A 'best' model from Whitehouse et al. (2012). Suitable
%           only for Antarctica.
%       The input can also be the path to the model file.
%		The default model is 'Steffen_ice6g_vm5a'.
%   r - The angular extent of the spherical cap radius in degrees
%   domain - The domain of interest
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A cell array of the domain of interest (string) and the buffer
%           value (scalar): {domain, buf}.
%       - An N-by-2 matrix of the domain of interest in [lon, lat] format.
%       - A GeoDomain object.
%   L - The bandwith of the spherical harmonic coefficients
%       - A scalar L, denoting the maximum degree of the SH coefficients.
%       - A 1-by-2 vector [Lmin, Lmax], denoting the minimum and maximum
%           degree of the SH coefficients.
%   days - Number of days to calculate the GIA change for
%   time - Vector of dates to calculate the GIA change for
%		The input can be in DATENUM or DATETIME format.
%   L - Maximum degree of the GIA model
%		If empty, the model is not truncated.
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the centre of the tapers in degrees
%       The default values are 0.
%   BeQuiet - Whether to surpress output messages
%		The default option is false.
%
% Output arguments
%   time - Time vector
%   slept - GIA correction time series projected onto the Slepian basis
%   sleptU, sleptL - Upper and lower bounds of the GIA correction time
%       series, if available
%   total - Total GIA correction
%
% Last modified by
%   2024/08/20, williameclee@arizona.edu (@williameclee)

function varargout = gia2slept(varargin)
    %% Initialisation
    % Parse inputs
    [time, model, domain, L, phi, theta, omega, beQuiet] = ...
        parseinputs(varargin{:});

    %% Loading the model
    % Get the yearly trend
    if strcmp(model, 'mascon')
        slept = mascon2slept('gia', domain, L, ...
            [time(1), time(end)] + [-1, 1]);
        [G, ~, ~, ~, N] = glmalpha_new(domain, L, "BeQuiet", beQuiet);
        hasBounds = false;
    else

        warning('off', 'SLEPIAN:gia2plmt:noBoundsToReturn');

        [plm, plmU, plmL] = gia2plmt( ...
            [], model, L, "BeQuiet", beQuiet);
        hasBounds = ~isnan(plmU) && ~isnan(plmL);

        % plm(:, 3:4) = plm(:, 3:4) / sqrt(4 * pi); % TEST
        L = conddefval(L, max(plm(:, 1)));

        %% Computing the basis
        if isa(domain, 'GeoDomain') || ismatrix(domain)
            [falpha, ~, N, ~, G] = plm2slep_new( ...
                plm, domain, L, "BeQuiet", beQuiet);

            if hasBounds
                falphaU = plm2slep_new( ...
                    plmU, domain, L, "BeQuiet", beQuiet);
                falphaL = plm2slep_new( ...
                    plmL, domain, L, "BeQuiet", beQuiet);
            end

        else
            [falpha, ~, N, ~, G] = plm2slep_new( ...
                plm, domain, L, phi, theta, omega, "BeQuiet", beQuiet);

            if hasBounds
                falphaU = plm2slep_new(plmU, domain, L, ...
                    phi, theta, omega, "BeQuiet", beQuiet);
                falphaL = plm2slep_new(plmL, domain, L, ...
                    phi, theta, omega, "BeQuiet", beQuiet);
            end

        end

        %% Getting the trend
        if isempty(time) || isscalar(time)
            % If time is scalar, interpret it as the day change
            if isempty(time)
                time = 365;
                deltaYear = 1;
            else
                deltaYear = time / 365;
            end

        else
            deltaYear = (time - time(1)) / 365;
        end

        slept = deltaYear(:) * falpha(:)';

        if hasBounds
            sleptU = deltaYear(:) * falphaU(:)';
            sleptL = deltaYear(:) * falphaL(:)';
        end

    end

    %% Getting the total
    truncation = round(N);

    if isa(domain, 'GeoDomain') || ismatrix(domain)
        eigfunINT = integratebasis_new( ...
            G, domain, truncation, "BeQuiet", beQuiet);
    else
        eigfunINT = integratebasis_new( ...
            G, domain, truncation, phi, theta, "BeQuiet", beQuiet);
    end

    eigfunINT = eigfunINT * (4 * pi * 6370e3 ^ 2) / 1e12;
    total = slept(:, 1:truncation) * eigfunINT';

    %% Collecting outputs and plotting
    if ~hasBounds
        sleptU = nan;
        sleptL = nan;
    end

    varargout = {time, slept, sleptU, sleptL, total};

    if nargout > 0
        return
    end

    plotgiamap(model, plm, falpha, deltaYear, domain, L);

end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Fallback values
    modelD = 'Steffen_ice6g_vm5a';
    domainD = {'greenland', 0.5};
    phiD = 0;
    thetaD = 0;
    omegaD = 0;

    % Allow skipping the time argument
    if nargin > 0 && ...
            (ischar(varargin{1}) || isstring(varargin{1}) || ...
            iscell(varargin{1}))
        varargin(2:end + 1) = varargin;
        varargin{1} = [];
    end

    % Parsing inputs
    p = inputParser;
    addOptional(p, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'Domain', domainD, ...
        @(x) ischar(x) || isstring(x) || iscell(x) ...
        || isa(x, 'GeoDomain') || isnumeric(x) || isempty(x));
    addOptional(p, 'L', [], ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'phi', phiD, @(x) isnumeric);
    addOptional(p, 'theta', thetaD, @(x) isnumeric);
    addOptional(p, 'omega', omegaD, @(x) isnumeric);
    addParameter(p, 'BeQuiet', 0.5, @(x) islogical(x) || isnumeric(x));

    parse(p, varargin{:});
    time = p.Results.Time(:);
    model = conddefval(p.Results.Model, modelD);
    domain = conddefval(p.Results.Domain, domainD);
    L = p.Results.L;
    phi = conddefval(p.Results.phi, phiD);
    theta = conddefval(p.Results.theta, thetaD);
    omega = conddefval(p.Results.omega, omegaD);
    beQuiet = uint8(double(p.Results.BeQuiet) * 2);

    if isdatetime(time)
        time = datenum(time); %#ok<DATNM>
    end

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    varargout = {time, model, domain, L, phi, theta, omega, beQuiet};

end

function plotgiamap(model, plm, slept, deltaYear, domain, L)
    %% Preparing the data
    meshSize = 1;
    mesh = plm2xyz(plm, meshSize, "BeQuiet", true);
    plmLcl = slep2plm_new(slept, domain, L, "BeQuiet", true);
    [meshLcl, lon, lat] = plm2xyz(plmLcl, meshSize, "BeQuiet", true);

    deltaYear = deltaYear(end);
    mesh = mesh * deltaYear;
    meshLcl = meshLcl * deltaYear;

    coastLonlat = gshhscoastline('c', 'LonOrigin', 180, "BeQuiet", true);
    domainLonlat = domain.Lonlat('LonOrigin', 180);

    [cLim, cStep] = optimalclim(meshLcl, 'Percentile', 1);
    mesh = max(min(mesh, cLim(2)), cLim(1));
    meshLcl = max(min(meshLcl, cLim(2)), cLim(1));

    %% Plotting
    % Protect underscore in model name
    if iscell(model)
        model = strjoin(model, '_');
    end

    model = strrep(model, '_', '\_');
    figure(999)
    set(gcf, "NumberTitle", 'off', "Name", ...
        sprintf('Localised GIA change in %.1f year(s) (%s)', ...
        deltaYear, upper(mfilename)))
    clf

    subplot(1, 2, 1)
    title(sprintf('Model: %s (global)', model))
    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", 'GIA change [kg/m^2]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, mesh, cLevels, 'LineStyle', 'none')
    plotqdm(coastLonlat, 'k');
    plotqdm(domainLonlat, 'k', 'LineWidth', 1);
    hold off

    subplot(1, 2, 2)
    title(sprintf('Model: %s (localised)', model))
    loadcbar(cLim, cStep, ...
        "Title", 'GIA change [kg/m^2]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, meshLcl, cLevels, 'LineStyle', 'none')
    plotqdm(coastLonlat, 'k');
    plotqdm(domainLonlat, 'k', 'LineWidth', 1);
    hold off
end
