%% PLM2SLEP
% Finds the spherical harmonic expansion coefficients into a SINGLE-CAP,
% potentially rotated, Slepian basis of a function whose real spherical
% harmonic expansion coefficients are given.
%
% Syntax
%   plm2slep(demoId)
%       Runs a demo with the specified name.
%   [falpha, V, N] = plm2slep(Plm, r, L, phi, theta, omega)
%       Finds the expansion coefficients of the function into the Slepian
%       basis of a polar cap.
%   [falpha, V, N] = plm2slep(Plm, domain, L)
%       Finds the expansion coefficients of the function into the Slepian
%       basis of a geographic domain.
%   [falpha, V, N] = plm2slep(__, nosort, truncation)
%       Finds the expansion with the specified sorting and the number of
%       eigenfunctions to use.
%   [__, MTAP, G] = plm2slep(__)
%
% Input arguments
%   Plm - Standard-type real spherical harmonic expansion coefficients
%   r - Radius of the concentration region in degrees
%       The default value is 30 degrees
%   domain - Geographic domain or a latitude-longitude pair
%       - A geographic domain (GeoDomain object).
%       - A string of the domain name.
%       - A cell array of the form {'domain name', buf}.
%       - A N-by-2 matrix of longitude-latitude vertices.
%   L - Bandwidth of the window
%       The default value is 18.
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the centre of the tapers in degrees
%       The default values are 0.
%   nosort - Whether to sort the eigenvalues and eigenfunctions
%       - false: Will sort the output according to the global eigenvalue.
%       - true: Will not sort thus the "block" sorting is preserved,
%       The default value is false (i.e. sort).
%   truncation - Number of largest eigenfunctions in which to expand
%       The default value is all the (L+1)^2 eigenfunctions.
%
% Output arguments
%   falpha - The expansion coefficients of the function into the Slepian
%       basis
%   V - The eigenvalues of the Slepian functions in question
%   N - The Shannon number
%   MTAP - The orders of the Slepian functions in question, if preserved
%   G - The matrix with the spherical harmonic expansion
%            coefficients of the Slepian functions used
%
% Examples
%   Check that single-order functions correctly transform back
%   >>  plm2slep('demo1')
%   Slight variation on the above with multiple same-order ones
%   >>  plm2slep('demo2')
%   Moon without the South Pole Aitken Basin
%   >>  plm2slep('demo3')
%
% See also
%   PTOSLEP, GLMALPHA, GLMALPHAPTO, SLEP2PLM
%
% Last modified by
%   2024/08/13, williameclee@arizona.edu (@williameclee)
%   2024/07/12, williameclee@arizona.edu (@williameclee)
%   2023/09/26, fjsimons@alum.mit.edu (@fjsimons)
%   2013/04/24, charig@princeton.edu (@charig)

function varargout = plm2slep_new(varargin)
    %% Initialisation
    % Add path to the auxiliary functions
    addpath(fullfile(fileparts(mfilename('fullpath')), 'demos'));
    % Demos
    if ischar(varargin{1}) || isstring(varargin{1})
        demoId = varargin{1};

        switch demoId
            case 'demo1'
                plm2slep_demo1
            case 'demo2'
                plm2slep_demo2
            case 'demo3'
                plm2slep_demo3
            otherwise
                error('Unknown demo ''%s''', upper(demoId))
        end

        return
    end

    % Parse inputs
    [lmcosi, domain, L, phi, theta, omega, nosort, J, beQuiet, GVN] = ...
        parseinputs(varargin{:});

    maxL = max(L);
    % The spherical harmonic dimension
    if isscalar(L)
        ldim = (L + 1) ^ 2;
    else
        ldim = (L(2) + 1) ^ 2 - L(1) ^ 2;
    end

    J = conddefval(J, ldim);

    if lmcosi(1) ~= 0 || lmcosi(2) ~= 1
        error('Spherical harmonics must start from degree zero')
    end

    %% Computing the projection
    % If it is the standard North-Polar cap or a geographic region, it's easy
    if ~isempty(GVN)
        G = GVN{1};
        V = GVN{2};
        N = GVN{3};
        MTAP = nan;
    elseif phi == 0 && theta == 0 && omega == 0
        % Get the Slepian basis; definitely not block-sorted as for the rotated
        % versions this will make no sense at all anymore
        % Glmalpha can handle a string, cell, or coordinates as TH, so this is ok
        [G, V, ~, ~, N, ~, MTAP, ~] = ...
            glmalpha_new(domain, L, "J", J, "BeQuiet", beQuiet);

    else
        % Need to get a complete GLMALPHA but for the rotated basis
        % Definitely, "single-order" has lost its meaning here, but the MTAP
        % will still identify what the order of the unrotated original was
        [G, V, ~, ~, N, ~, MTAP, ~] = ...
            glmalphapto(domain, L, phi, theta, omega);
    end

    if ~nosort
        % Sort by decreasing eigenvalue
        [V, vi] = sort(V, 'descend');
        G = G(:, vi);

        if ~isnan(MTAP)
            MTAP = MTAP(vi);
        end

        % If you don't do this, the eigenfunctions are ordered in the way
        % that they correspond to single-orders back when, unrotated, they
        % belonged to a polar cap, and the eigenvalues are sorted within
        % these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.
    end

    % Get the mapping from LMCOSI into not-block-sorted GLMALPHA
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ronm] = addmon(maxL);

    % Make sure that the requested L acts as truncation on lmcosi
    % or if we don't have enough, pad with zeros
    if size(lmcosi, 1) < addmup(maxL)
        [~, ~, ~, lmcosipad] = addmon(maxL);
        lmcosi = [lmcosi; lmcosipad(size(lmcosi, 1) + 1:end, :)];
    else
        lmcosi = lmcosi(1:addmup(maxL), :);
    end

    % Perform the expansion of the signal into the Slepian basis
    falpha = G' * lmcosi(2 * size(lmcosi, 1) + ronm(1:(maxL + 1) ^ 2));

    %% Returning requested outputs
    varargout = {falpha, V, N, MTAP, G};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    domainD = 30;
    LD = 18;
    phiD = 0;
    thetaD = 0;
    omegaD = 0;
    nosortD = false;
    JD = [];
    upscaleD = 0;

    p = inputParser;
    addRequired(p, 'lmcosi', ...
        @(x) isnumeric(x) || ischar(x));
    addOptional(p, 'Domain', domainD, ...
        @(x) ischar(x) || iscell(x) || isa(x, "GeoDomain") || ...
        isnumeric(x) || isempty(x));
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'phi', phiD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'theta', thetaD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'omega', omegaD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'nosort', nosortD, ...
        @(x) isnumeric(x) || islogical(x) || isempty(x));
    addOptional(p, 'J', JD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Upscale', upscaleD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'MoreDomainSpecs', {}, @iscell);
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'GVN', {}, @(x) iscell(x) && length(x) == 3);
    parse(p, varargin{:});

    lmcosi = p.Results.lmcosi;
    domain = conddefval(p.Results.Domain, domainD);
    L = conddefval(p.Results.L, LD);
    phi = conddefval(p.Results.phi, phiD);
    theta = conddefval(p.Results.theta, thetaD);
    omega = conddefval(p.Results.omega, omegaD);
    nosort = conddefval(p.Results.nosort, nosortD);
    J = conddefval(p.Results.J, JD);
    upscale = conddefval(p.Results.Upscale, upscaleD);
    moreRegionSpecs = p.Results.MoreDomainSpecs;
    beQuiet = logical(p.Results.BeQuiet);
    GVN = p.Results.GVN;

    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, "Upscale", upscale, moreRegionSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = GeoDomain(domain{1}, ...
            "Buffer", domain{2}, "Upscale", upscale, moreRegionSpecs{:});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{1}, "Upscale", upscale, domain{2:end});
    end

    varargout = {lmcosi, domain, L, phi, theta, omega, nosort, J, beQuiet, GVN};
end
