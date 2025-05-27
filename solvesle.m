%% SOLVESLE
% Solves the (elastic) sea level equation (SLE) for a given forcing and
% ocean domain.
% The ocean function is assumed constant.
%
% Syntax
%   rslPlm = solvesle(forcing)
%   rslPlm = solvesle(forcing, forcingStd, L)
%   rslPlm = solvesle(__, "Name", Value)
%   [rslPlm, rslHatPlm, barystatic] = solvesle(__)
%
% Input arguments
%   forcing - Forcing spherical harmonic coefficients
%       The first two columns (degree and order) can be omitted, assuming the data is in the form lmcosi.
%       The forcing will not be localised, i.e. signal in the ocean is retained even when enforcing mass conservation (could be useful for, e.g. GAD).
%   forcingStd - Standard deviation of the forcing coefficients
%   L - Bandwidth of the data
%       The default degree is the degree of forcing.
%   Ocean - Ocean domain
%       A geographic domain (GeoDomain object).
%       The default domain is all oceans with a buffer of 0.5°.
%   Frame - Reference frame
%       Centre of mass (CM) frame or centre of figure (CF) frame.
%       The default frame is the CF frame, which is the frame used in slepian_delta
%   RotationFeedback - Whether to include rotation feedback
%       The default value is true
%   MaxIter - Maximum number of iterations
%       The default value is 10
%   OceanKernel - Precomputed ocean kernel
%
% Output arguments
%   rslLoadPlm - Mass load spherical harmonic coefficients
%       Unit: kg/m^2 (<=> mm fresh water)
%   rslLoadLclPlm - Localised sea level spherical harmonic coefficients
%       Unit: kg/m^2 (<=> mm fresh water)
%   barystatic - Static ocean mass change
%       Unit: mm
%
% Notes
%   In the implementation, all SH coefficient variables with a 's' suffix
%   means they are flattened to a column vector.
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/05/27, williameclee@arizona.edu (@williameclee)

function [rslLoadPlm, rslOceanLoadPlm, barystaticSL] = solvesle(varargin)
    %% Initialisation
    % Parse input arguments
    [forcingLoadPlm, forcingLoadStdPlm, includesDO, computeError, ...
         L, ocean, frame, doRotationFeedback, maxIter, beQuiet, ...
         oceanKernel, oceanFunPlm, kernelOrder] = parseinputs(varargin{:});

    % Truncate the input to the correct degree
    if size(forcingLoadPlm, 1) < addmup(L)
        forcingLoadPlm(addmup(L), 2) = 0;
        forcingLoadStdPlm(addmup(L), 2) = 0;
    elseif size(forcingLoadPlm, 1) > addmup(L)
        forcingLoadPlm = forcingLoadPlm(1:addmup(L), :);
        forcingLoadStdPlm = forcingLoadStdPlm(1:addmup(L), :);
    end

    if isempty(oceanFunPlm)

        if ~beQuiet
            counterTxt = sprintf('SLE solver progress: Computing ocean function %3d%% (%2d/%2d iterations)', ...
                0, 0, maxIter);
            fprintf('%s\n', counterTxt)
            counterTxtLen = length(counterTxt);
        end

        [~, ~, ~, ~, ~, oceanFunPlm] = geoboxcap(L, ocean, "BeQuiet", true);

        if ~beQuiet
            counterTxt = sprintf('SLE solver progress: Computed ocean function  %3d%% (%2d/%2d iterations)', ...
                0, 0, maxIter);
            fprintf('%s%s\n', sprintf(repmat('\b', 1, counterTxtLen + 1)), counterTxt);
            counterTxtLen = length(counterTxt);
        end

    else
        counterTxtLen = -1;
    end

    if size(oceanFunPlm, 2) == 4
        oceanFunPlm = oceanFunPlm(:, 3:4);
    end

    if isempty(oceanKernel)
        oceanKernel = kernelcp_new(L, ocean, "BeQuiet", true);
    end

    if isempty(kernelOrder)
        kernelOrder = kernelorder(L);
    end

    kernelOrder = kernelOrder(:);

    forcingLoadPlms = forcingLoadPlm(kernelOrder);
    forcingLoadStdPlms = forcingLoadStdPlm(kernelOrder);
    oceanFunPlms = oceanFunPlm(kernelOrder);

    %% Constants
    % Many values from Adhikari et al. (2016; doi:10.5194/gmd-9-1087-2016)
    WATER_DENSITY = 1000; % to be consistent with other functions, mainly MASS2WEQ
    EARTH_DENSITY = 5517; % to be consistent with PLM2POT
    EARTH_RADIUS = 6371e3;
    EARTH_ANGULAR_VELOCITY = 7.2921e-5;
    EQUATORIAL_INERTIA = 8.0077e37;
    POLAR_INERTIA = 8.0345e37;
    CHANDLER_WOBBLE_FREQUENCY = 2.4405e-7;
    GRAVITY = 9.81;

    [~, degree] = addmon(L);
    degrees = [degree, degree];
    degrees = degrees(kernelOrder);
    llnGeoids = lovenumber(degrees, 'loadinggravitationalpotential', frame);
    llnVlms = lovenumber(degrees, 'loadingverticaldisplacement', frame);
    llnFactors = 1 + llnGeoids - llnVlms;
    tlnGeoids = lovenumber(degrees, 'tidalgravitationalpotential', frame);
    tlnVlms = lovenumber(degrees, 'tidalverticaldisplacement', frame);
    tlnFactors = 1 + tlnGeoids - tlnVlms;
    llnGeoid2 = lovenumber(2, 'loadinggravitationalpotential', frame);

    %% Iteration
    barystaticSL =- forcingLoadPlms(1, :) / oceanFunPlms(1, :) / WATER_DENSITY;
    rslPlms = oceanFunPlms * barystaticSL;
    rslOceanPlms = rslPlms;

    for iIter = 1:maxIter
        loadPlms = forcingLoadPlms + rslOceanPlms * WATER_DENSITY;

        rslGravityPlms = 3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1) .* loadPlms;

        if doRotationFeedback
            rotPotPlms = zeros(size(rslGravityPlms));
            % C00
            rotPotPlms(1, :) = 2/3 * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 ...
                * (- (1 + llnGeoid2) / POLAR_INERTIA) * 8 * pi / 3 * EARTH_RADIUS ^ 4 ...
                * (loadPlms(1, :) - loadPlms(5, :) / sqrt(5));
            % C20
            rotPotPlms(5, :) = -rotPotPlms(1, :) / sqrt(5);
            % C21/S21
            rotPotPlms(6:7, :) = -1 / sqrt(15) * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 ...
                * EARTH_ANGULAR_VELOCITY * (1 + llnGeoid2) / (EQUATORIAL_INERTIA * CHANDLER_WOBBLE_FREQUENCY) ...
                * (-4 * pi / sqrt(15)) * EARTH_RADIUS ^ 4 * loadPlms(6:7, :);
            rslRotPlms = tlnFactors / GRAVITY .* rotPotPlms;
            rslPlms = rslGravityPlms + rslRotPlms;
        else
            rslPlms = rslGravityPlms;
        end

        rslOceanPlms = oceanKernel * rslPlms;

        % Enforce mass conservation
        rslOceanPlms(1, :) = -forcingLoadPlms(1, :) / WATER_DENSITY;
        rslPlms(1, :) = rslOceanPlms(1, :) / oceanFunPlms(1, :);

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
    rslLoadPlms = rslPlms * WATER_DENSITY;
    rslOceanLoadPlms = rslOceanPlms * WATER_DENSITY;

    barystaticSL = rslOceanPlms(1, :) ./ oceanFunPlms(1, :) * 1e3; % m -> mm

    rslLoadPlm = zeros([addmup(L), 2]);
    rslOceanLoadPlm = zeros([addmup(L), 2]);
    rslLoadPlm(kernelOrder) = rslLoadPlms;
    rslOceanLoadPlm(kernelOrder) = rslOceanLoadPlms;

    if includesDO
        [order, degree] = addmon(L);
        rslLoadPlm = [degree, order, rslLoadPlm];
        rslOceanLoadPlm = [degree, order, rslOceanLoadPlm];
    end

    % Clean counter
    if ~beQuiet
        fprintf('%s', sprintf(repmat('\b', 1, counterTxtLen + 1)));
    end

end

%% Subfunctions
% Parse input arguments
function varargout = parseinputs(varargin)
    ip = inputParser;
    addRequired(ip, 'ForcingPlm', ...
        @(x) isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2));
    addOptional(ip, 'ForcingStdPlm', [], ...
        @(x) (isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2)) || isempty(x));
    addOptional(ip, 'L', [], ...
        @(x) (isnumeric(x) && isscalar(x) && x > 0) || isempty(x));
    addOptional(ip, 'Ocean', GeoDomain('alloceans', "Buffer", 0.5), ...
        @(x) isa(x, 'GeoDomain') || (isnumeric(x) && ismatrix(x) && size(x, 2) == 2));
    addOptional(ip, 'frame', 'CF', ...
        @(x) ischar(validatestring(upper(x), {'CM', 'CF'})));
    addOptional(ip, 'RotationFeedback', true, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'maxIter', 10, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'OceanKernel', [], @(x) isnumeric(x));
    addParameter(ip, 'OceanFunction', [], @(x) isnumeric(x));
    addParameter(ip, 'KernelOrder', [], @(x) isnumeric(x));
    parse(ip, varargin{:});
    forcingPlm = squeeze(ip.Results.ForcingPlm);
    forcingStdPlm = squeeze(ip.Results.ForcingStdPlm);
    L = ip.Results.L;
    ocean = ip.Results.Ocean;
    frame = upper(ip.Results.frame);
    doRotationFeedback = ip.Results.RotationFeedback;
    maxIter = ip.Results.maxIter;
    beQuiet = ip.Results.BeQuiet;
    oceanKernel = ip.Results.OceanKernel;
    oceanFunSph = ip.Results.OceanFunction;
    kernelOrder = ip.Results.KernelOrder;

    % Check input sizes
    includesDO = size(forcingPlm, 2) == 4;

    if isempty(forcingStdPlm) || ...
            all(forcingStdPlm(:, end - 1:end) == 0, "all")
        computeError = false;
    else

        if ~isequal(size(forcingStdPlm), size(forcingPlm))
            error( ...
                sprintf('%s:InvalidInput:InputSizeNotMatch', upper(mfilename)), ...
                'Input size of FORCINGSTDPLM ([%s]) does not match FORCINGPLM ([%s])', ...
                strjoin(string(size(forcingStdPlm)), ', '), ...
                strjoin(string(size(forcingPlm)), ', '));
        end

        computeError = true;

    end

    % Find the degree of the spherical harmonic coefficients
    if includesDO

        if isempty(L)
            L = max(forcingPlm(:, 1));
        end

        forcingPlm = forcingPlm(:, 3:4);
    elseif isempty(L)
        L = finddegree(forcingPlm);
    end

    varargout = {forcingPlm, forcingStdPlm, includesDO, computeError, L, ocean, frame, doRotationFeedback, maxIter, beQuiet, oceanKernel, oceanFunSph, kernelOrder};
end

% Find the degree of the spherical harmonic coefficients based on the number of terms
function L = finddegree(Plm)
    terms = size(Plm, 1);
    syms x y
    eq = x ^ 2 + 3 * x + (2 - 2 * y) == 0;
    eq = subs(eq, y, terms);
    Ls = solve(eq, x);
    L = Ls(Ls > 0);
end
