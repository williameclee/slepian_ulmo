%% SOLVESLE
% Solves the sea level equation (SLE) for a given forcing and ocean domain.
% Note this function only solves the elastic component, and assumes the ocean function is constant.
%
% Syntax
%   rslPlm = solvesle(forcingPlm)
%   rslPlm = solvesle(forcingPlm, L)
%   rslPlm = solvesle(__, "Name", Value)
%   [rslPlm, rslHatPlm, barystatic] = solvesle(__)
%
% Input arguments
%   forcingPlm - Forcing spherical harmonic coefficients
%       The first two columns (degree and order) can be omitted, assuming the data is in the form lmcosi.
%   L - Bandwidth of the window
%       The default degree is 96
%   "Ocean" - Ocean domain
%       A geographic domain (GeoDomain object).
%       The default domain is all oceans with a buffer of 1 degree.
%   "Frame" - Reference frame
%       Centre of mass (CM) frame or centre of figure (CF) frame.
%       The default frame is the CF frame, which is the frame used in slepian_delta
%   "RotationFeedback" - Whether to include rotation feedback
%       The default value is true
%   "MaxIter" - Maximum number of iterations
%       The default value is 10
%   "OceanKernel" - Precomputed ocean kernel
%
% Output arguments
%   rslPlm - Sea level spherical harmonic coefficients
%   rslHatPlm - Localised sea level spherical harmonic coefficients
%   barystatic - Static ocean mass change
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)

function [rslPlm, rslHatPlm, barystatic] = solvesle(forcingPlm, varargin)
    %% Initialisation
    ip = inputParser;
    addRequired(ip, 'forcingPlm', ...
        @(x) isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2));
    % addOptional(ip, 'giaPlm', ...
    %     @(x) (isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2)) || isempty(x));
    addOptional(ip, 'L', [], ...
        @(x) (isnumeric(x) && isscalar(x) && x > 0) || isempty(x));
    addOptional(ip, 'ocean', GeoDomain('alloceans'), ...
        @(x) isa(x, 'GeoDomain') || (isnumeric(x) && ismatrix(x) && size(x, 2) == 2));
    addOptional(ip, 'frame', 'CF', ...
        @(x) ischar(validatestring(upper(x), {'CM', 'CF'})));
    addOptional(ip, 'RotationFeedback', true, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'maxIter', 10, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'OceanKernel', [], @(x) isnumeric(x));
    parse(ip, forcingPlm, varargin{:});
    forcingPlm = ip.Results.forcingPlm;
    % giaPlm = ip.Results.giaPlm;
    L = ip.Results.L;
    ocean = ip.Results.ocean;
    frame = upper(ip.Results.frame);
    doRotationFeedback = ip.Results.RotationFeedback;
    maxIter = ip.Results.maxIter;
    beQuiet = ip.Results.BeQuiet;
    Kocean = ip.Results.OceanKernel;

    % if ~isempty(giaPlm)
    %     doGia = true;

    %     if ~isequal(size(forcingPlm), size(giaPlm))
    %         error('Forcing and GIA must have the same size')
    %     end

    % else
    %     doGia = false;
    % end

    includesDO = size(forcingPlm, 2) == 4;

    if includesDO

        if isempty(L)
            L = max(forcingPlm(:, 1));
        end

        forcingPlm = forcingPlm(:, 3:4);

        % if doGia
        %     giaPlm = giaPlm(:, 3:4);
        % end

    elseif isempty(L)
        L = finddegree(forcingPlm);
    end

    if size(forcingPlm, 1) < addmup(L)
        forcingPlm(addmup(L), 2) = 0;
        % giaPlm(addmup(L), 2) = 0;
    elseif size(forcingPlm, 1) > addmup(L)
        forcingPlm = forcingPlm(1:addmup(L), :);
        % giaPlm = giaPlm(1:addmup(L), :);
    end

    % forcingPlm = localise(forcingPlm, ocean, L, "Inverse", true);
    % giaLandPlm = localise(giaPlm, ocean, L, "Inverse", true);
    % forcingPlm = forcingPlm - giaLandPlm;

    if ~beQuiet
        counterTxt = sprintf('SLE solver progress: Computing ocean function %3d%% (%2d/%2d iterations)', ...
            0, 0, maxIter);
        fprintf('%s\n', counterTxt)
        counterTxtLen = length(counterTxt);
    end

    [~, ~, ~, ~, ~, oceanPlm] = geoboxcap(L, ocean, "BeQuiet", true);
    oceanPlm = oceanPlm(:, 3:4);

    if ~beQuiet
        counterTxt = sprintf('SLE solver progress: Computed ocean function  %3d%% (%2d/%2d iterations)', ...
            0, 0, maxIter);
        fprintf('%s%s\n', sprintf(repmat('\b', 1, counterTxtLen + 1)), counterTxt);
        counterTxtLen = length(counterTxt);
    end

    %% Constants
    % Many values from Adhikari et al. (2016; doi:10.5194/gmd-9-1087-2016)
    WATER_DENSITY = 1028;
    EARTH_DENSITY = 5515;
    EARTH_RADIUS = 6371e3;
    EARTH_ANGULAR_VELOCITY = 7.2921e-5;
    EQUATORIAL_INERTIA = 8.0077e37;
    POLAR_INERTIA = 8.0345e37;
    CHANDLER_WOBBLE_FREQUENCY = 2.4405e-7;
    GRAVITATIONAL_ACCELERATION = 9.81;

    [~, degree] = addmon(L);
    llngeoid = lovenumber(degree, 'loadinggravitationalpotential', frame);
    llnvdisp = lovenumber(degree, 'loadingverticaldisplacement', frame);
    llnfactor = 1 + llngeoid - llnvdisp;
    tlngeoid = lovenumber(degree, 'tidalgravitationalpotential', frame);
    tlnvdisp = lovenumber(degree, 'tidalverticaldisplacement', frame);
    tlnfactor = 1 + tlngeoid - tlnvdisp;

    %% Iteration
    barystatic =- forcingPlm(1, 1) / oceanPlm(1, 1);

    rslPlm = oceanPlm * barystatic; % is actually density anomaly
    rslHatPlm = rslPlm;

    if isempty(Kocean)
        Kocean = kernelcp_new(L, ocean, "BeQuiet", true);
    end

    for iIter = 1:maxIter
        loadPlm = forcingPlm + rslHatPlm;

        % if doGia
        %     loadPlm = loadPlm + giaPlm;
        % end

        rslLoadPlm = 3 * WATER_DENSITY / EARTH_DENSITY * llnfactor ./ (2 * degree + 1) .* loadPlm;

        if doRotationFeedback
            rotPotPlm = zeros([addmup(L), 2]);
            rotPotPlm(1, 1) = 2/3 * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 * (- (1 + lovenumber(2, 'loadinggravitationalpotential', frame)) / POLAR_INERTIA) * 8 * pi / 3 * EARTH_RADIUS ^ 4 * (loadPlm(1, 1) - loadPlm(4, 1) / sqrt(5));
            rotPotPlm(4, 1) = -rotPotPlm(1, 1) / sqrt(5);
            rotPotPlm(5, :) = -1 / sqrt(15) * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 * EARTH_ANGULAR_VELOCITY * (1 + lovenumber(2, 'loadinggravitationalpotential', frame)) / (EQUATORIAL_INERTIA * CHANDLER_WOBBLE_FREQUENCY) * (-4 * pi / sqrt(15)) * EARTH_RADIUS ^ 4 * loadPlm(5, :);
            rslRotPlm = WATER_DENSITY * tlnfactor / GRAVITATIONAL_ACCELERATION .* rotPotPlm;
            rslPlm = rslLoadPlm + rslRotPlm;
        else
            rslPlm = rslLoadPlm;
        end

        heightDiff =- (forcingPlm(1, 1) / oceanPlm(1, 1)) - ...
            sum(rslPlm .* oceanPlm / oceanPlm(1, 1), 'all');
        % still actually in density anomaly
        rslPlm(1, 1) = rslPlm(1, 1) + heightDiff;
        rslHatPlm = localise(rslPlm, ocean, L, "K", Kocean);

        % Counter
        if ~beQuiet
            div = 24;
            progress = floor(iIter / maxIter * div);
            counterTxt = sprintf('SLE solver progress: %s%s %3d%% (%2d/%2d iterations)', ...
                repmat('=', 1, progress), repmat(' ', 1, div - progress), iIter / maxIter * 100, iIter, maxIter);
            fprintf('%s%s\n', sprintf(repmat('\b', 1, counterTxtLen + 1)), counterTxt);
            counterTxtLen = length(counterTxt);
        end

    end

    %% Post-processing
    barystatic = rslHatPlm(1, 1) / oceanPlm(1, 1);

    if includesDO
        [order, degree] = addmon(L);
        rslPlm = [degree, order, rslPlm];
        rslHatPlm = [degree, order, rslHatPlm];
    end

end

function L = finddegree(Plm)
    terms = size(Plm, 1);
    syms x y
    eq = x ^ 2 + 3 * x + (2 - 2 * y) == 0;
    eq = subs(eq, y, terms);
    Ls = solve(eq, x);
    L = Ls(Ls > 0);
end
