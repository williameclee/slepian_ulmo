%% SOLVESLE
% Solves the sea level equation (SLE) for a given forcing and ocean domain.
% Note this function only solves the elastic component of the modern load, and the viscous response is from a GIA model (if specified).
% The ocean function is assumed constant.
%
% Syntax
%   rslPlm = solvesle(forcingSph)
%   rslPlm = solvesle(forcingSph, L)
%   rslPlm = solvesle(__, "Name", Value)
%   [rslPlm, rslHatPlm, barystatic] = solvesle(__)
%
% Input arguments
%   forcingSph - Forcing spherical harmonic coefficients
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
%   "GiaGeoidSph" & "GiaVlmSph" - GIA geoid and bedrock deformation spherical harmonic coefficients
%       If not specified, the GIA feedback is not included (which should be able to be added back later on without problem).
%   "StericSph" - Steric sea level spherical harmonic coefficients
%       If not specified, the steric component is not considered when converting between the load and the relative sea level.
%       Unit: mm
%   "MaxIter" - Maximum number of iterations
%       The default value is 10
%   "OceanKernel" - Precomputed ocean kernel
%
% Output arguments
%   rslLoadPlm - Mass load spherical harmonic coefficients
%       Unit: kg/m^2 (<=> mm fresh water)
%   rslLoadLclPlm - Localised sea level spherical harmonic coefficients
%       Unit: kg/m^2 (<=> mm fresh water)
%   barystatic - Static ocean mass change
%       Unit: mm
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/03/28, williameclee@arizona.edu (@williameclee)

function [rslLoadSph, rslLoadLclSph, barystaticSL] = solvesle(forcingSph, varargin)
    %% Initialisation
    ip = inputParser;
    addRequired(ip, 'ForcingSph', ...
        @(x) isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2));
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
    addOptional(ip, 'GiaGeoidSph', [], ...
        @(x) (isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2)) || isempty(x));
    addOptional(ip, 'GiaVlmSph', [], ...
        @(x) (isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2)) || isempty(x));
    % addOptional(ip, 'StericSph', [], ...
    %     @(x) (isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2)) || isempty(x));
    addOptional(ip, 'OtherForcingSph', [], ...
        @(x) (isnumeric(x) && ismatrix(x) && (size(x, 2) == 4 || size(x, 2) == 2)) || isempty(x));
    addParameter(ip, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'OceanKernel', [], @(x) isnumeric(x));
    addParameter(ip, 'OceanFunction', [], @(x) isnumeric(x));
    addParameter(ip, 'KernelOrder', [], @(x) isnumeric(x));
    parse(ip, forcingSph, varargin{:});
    forcingSph = squeeze(ip.Results.ForcingSph);
    % giaGeoidSph = squeeze(ip.Results.GiaGeoidSph);
    % giaVlmSph = squeeze(ip.Results.GiaVlmSph);
    % stericSph = squeeze(ip.Results.StericSph);
    otherForcingSph = squeeze(ip.Results.OtherForcingSph);
    L = ip.Results.L;
    ocean = ip.Results.Ocean;
    frame = upper(ip.Results.frame);
    doRotationFeedback = ip.Results.RotationFeedback;
    maxIter = ip.Results.maxIter;
    beQuiet = ip.Results.BeQuiet;
    oceanKernel = ip.Results.OceanKernel;
    oceanFunSph = ip.Results.OceanFunction;
    kernelOrder = ip.Results.KernelOrder;

    % if ~isempty(giaGeoidSph) && ~isempty(giaVlmSph)
    %     doGia = true;

    %     % Check if the forcing and GIA have the same size
    %     if ~isequal(size(forcingSph), size(giaGeoidSph))
    %         error('Forcing and GIA geoid deformation must have the same size')
    %     end

    %     if ~isequal(size(forcingSph), size(giaVlmSph))
    %         error('Forcing and GIA VLM must have the same size')
    %     end

    % else
    %     doGia = false;
    % end

    % if ~isempty(stericSph) && ...
    %         ~all(stericSph(:, end - 1:end) == 0, "all")
    %     doSteric = true;
    %     stericLclSph = localise(stericSph, ocean, L);
    % else
    %     doSteric = false;
    % end

    if ~isempty(otherForcingSph) && ...
            ~all(otherForcingSph(:, end - 1:end) == 0, "all")
        doOtherForcing = true;
    else
        doOtherForcing = false;
    end

    % Check if the forcing includes degree and order
    includesDO = size(forcingSph, 2) == 4;

    if includesDO

        if isempty(L)
            L = max(forcingSph(:, 1));
        end

        forcingSph = forcingSph(:, 3:4);

        % if doGia
        %     giaGeoidSph = giaGeoidSph(:, 3:4);
        %     giaVlmSph = giaVlmSph(:, 3:4);
        % end

        % if doSteric
        %     stericLclSph = stericLclSph(:, 3:4);
        % end

    elseif isempty(L)
        L = finddegree(forcingSph);
    end

    % if doSteric
    %     stericLclSph = stericLclSph / 1000; % mm -> m
    % end

    % Truncate the input to the correct degree
    if size(forcingSph, 1) < addmup(L)
        forcingSph(addmup(L), 2) = 0;
        %     giaGeoidSph(addmup(L), 2) = 0;
        %     giaVlmSph(addmup(L), 2) = 0;
    elseif size(forcingSph, 1) > addmup(L)
        forcingSph = forcingSph(1:addmup(L), :);
        % giaGeoidSph = giaGeoidSph(1:addmup(L), :);
    end

    if isempty(oceanFunSph)

        if ~beQuiet
            counterTxt = sprintf('SLE solver progress: Computing ocean function %3d%% (%2d/%2d iterations)', ...
                0, 0, maxIter);
            fprintf('%s\n', counterTxt)
            counterTxtLen = length(counterTxt);
        end

        [~, ~, ~, ~, ~, oceanFunSph] = geoboxcap(L, ocean, "BeQuiet", true);

        if ~beQuiet
            counterTxt = sprintf('SLE solver progress: Computed ocean function  %3d%% (%2d/%2d iterations)', ...
                0, 0, maxIter);
            fprintf('%s%s\n', sprintf(repmat('\b', 1, counterTxtLen + 1)), counterTxt);
            counterTxtLen = length(counterTxt);
        end

    else
        counterTxtLen = -1;
    end

    if size(oceanFunSph, 2) == 4
        oceanFunSph = oceanFunSph(:, 3:4);
    end

    if isempty(kernelOrder)
        kernelOrder = kernelorder(L);
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
    GRAVITY = 9.81;

    [~, degree] = addmon(L);
    llnGeoid = lovenumber(degree, 'loadinggravitationalpotential', frame);
    llnVlm = lovenumber(degree, 'loadingverticaldisplacement', frame);
    llnFactor = 1 + llnGeoid - llnVlm;
    tlnGeoid = lovenumber(degree, 'tidalgravitationalpotential', frame);
    tlnVlm = lovenumber(degree, 'tidalverticaldisplacement', frame);
    tlnFactor = 1 + tlnGeoid - tlnVlm;

    %% Iteration
    barystaticSL =- forcingSph(1, 1) / oceanFunSph(1, 1);

    rslSph = oceanFunSph * barystaticSL / WATER_DENSITY;
    rslLclSph = rslSph;

    if isempty(oceanKernel)
        oceanKernel = kernelcp_new(L, ocean, "BeQuiet", true);
    end

    for iIter = 1:maxIter

        % if doSteric
        %     loadSph = forcingSph + (rslLclSph - stericLclSph) * WATER_DENSITY;
        % else
        loadSph = forcingSph + rslLclSph * WATER_DENSITY;
        % end

        rslGravitySph = 3 / EARTH_DENSITY * llnFactor ./ (2 * degree + 1) .* loadSph;

        if doRotationFeedback
            rotPotSph = zeros([addmup(L), 2]);
            % C00
            rotPotSph(1, 1) = 2/3 * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 ...
                * (- (1 + llnGeoid(2)) / POLAR_INERTIA) * 8 * pi / 3 * EARTH_RADIUS ^ 4 * (loadSph(1, 1) - loadSph(4, 1) / sqrt(5));
            rotPotSph(4, 1) = -rotPotSph(1, 1) / sqrt(5);
            % C20
            rotPotSph(5, :) = -1 / sqrt(15) * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 ...
                * EARTH_ANGULAR_VELOCITY * (1 + llnGeoid(2)) / (EQUATORIAL_INERTIA * CHANDLER_WOBBLE_FREQUENCY) ...
                * (-4 * pi / sqrt(15)) * EARTH_RADIUS ^ 4 * loadSph(5, :);
            rslRotSph = tlnFactor / GRAVITY .* rotPotSph;
            rslSph = rslGravitySph + rslRotSph;
        else
            rslSph = rslGravitySph;
        end

        % if doGia
        %     rslSph = rslSph + (giaGeoidSph - giaVlmSph);
        % end

        if doOtherForcing
            rslSph = rslSph + otherForcingSph / WATER_DENSITY;
        end

        % % if doSteric
        % %     heightDiff = -(forcingSph(1, 1) / oceanFunSph(1, 1)) / WATER_DENSITY ...
        % %         -sum((rslSph - stericLclSph) .* oceanFunSph / oceanFunSph(1, 1), 'all');
        % % else
        % % end

        % heightDiff =- forcingSph(1, 1) / WATER_DENSITY / oceanFunSph(1, 1) ...
        %     -sum(rslSph .* oceanFunSph, 'all') / oceanFunSph(1, 1);
        % rslSph(1, 1) = rslSph(1, 1) + heightDiff;
        % rslLclSph = localise(rslSph, ocean, L, ...
        %     "K", oceanKernel, "KernelOrder", kernelOrder);
        rslLclSph(1, 1) = -forcingSph(1, 1) / WATER_DENSITY;
        rslSph(1, 1) = rslLclSph(1, 1) / oceanFunSph(1, 1);

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
    if doOtherForcing
        rslSph = rslSph - otherForcingSph / WATER_DENSITY;
        rslLclSph = localise(rslSph, ocean, L, "K", oceanKernel, "KernelOrder", kernelOrder);
    end

    % if doSteric
    %     rslLoadSph = (rslSph - stericLclSph) * WATER_DENSITY;
    %     rslLoadLclSph = (rslLclSph - stericLclSph) * WATER_DENSITY;
    % else
    rslLoadSph = rslSph * WATER_DENSITY;
    rslLoadLclSph = rslLclSph * WATER_DENSITY;
    % end

    barystaticSL = rslLclSph(1, 1) / oceanFunSph(1, 1) * 1000; % m -> mm

    if includesDO
        [order, degree] = addmon(L);
        rslLoadSph = [degree, order, rslLoadSph];
        rslLoadLclSph = [degree, order, rslLoadLclSph];
    end

    % Clean counter
    if ~beQuiet
        fprintf('%s', sprintf(repmat('\b', 1, counterTxtLen + 1)));
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
