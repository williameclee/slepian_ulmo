%% PLM2XYZ
% Inverse (4*pi-normalized real) spherical harmonic transform.
%
% Compute a spatial field from spherical harmonic coefficients given as
% [l m Ccos Csin] (not necessarily starting from zero, but sorted), with
% degree resolution 'h' [default: approximate Nyquist degree].
%
% Syntax
% [r, lon, lat, Plm, h] = plm2xyz(lmcosi, h)
% [r, lon, lat, Plm, h] = plm2xyz(lmcosi, h, c11cmn, lmax, latmax, Plm)
% [r, lon, lat, Plm, h] = plm2xyz(lmcosi, lat, lon, lmax, latmax, Plm)
%
% Input arguments
%   lmcosi - Matrix listing l,m,cosine and sine expansion coefficients
%       e.g. those coming out of ADDMON
%   h - Longitude/ latitude spacing, in degrees [default: Nyquist] OR
%       "lat": a column vector with latitudes [degrees]
%   c11cmn - Corner nodes of lon/lat grid [default: 0 90 360 -90] OR
%       "lon": a column vector with longitudes [degrees]
%   lmax - Maximum bandwidth expanded at a time
%       The default value is 720.
%   latmax - Maximum linear size of the latitude grid allowed
%       The default value is Inf.
%   Plm - The appropriate Legendre polynomials should you already have them
%
% Output arguments
%   r - The field (matrix for a grid, vector for scattered points)
%   lon, lat - The grid (matrix) or evaluation points (vector), in degrees
%   Plm - The set of appropriate Legendre polynomials should you want them
%   h - Longitude/latitude spacing, in degrees
%
% Examples
%   >>  plm2xyz('demo1') % Illustrates forward and inverse transform
%   >>  plm2xyz('demo2',fra) % with 'fra' a data fraction
%   >>  plm2xyz('demo3') % Plots EGM96 versus EGM2008 free-air anomaly
%   >>  plm2xyz('demo4') % Plots GTM3AR versus EGM2008 topography
%   >>  plm2xyz('demo5') % Plots EGM96 geopotential
%   >>  plm2xyz('demo6') % EGM2008 topography somewhere
%   >>  plm2xyz('demo7') % Tests the expansion to scattered points
%   >>  plm2xyz('demo8') % EGM2008 gravity globally
%
% Notes
%   The degree range is now split intelligently in blocks of degrees whose
%   memory requirements never exceed the initial expansion from 0 to lmax.
%   For very high bandwidth models, specify a small region, and this
%   together with 'h' will determine memory requirements. Built-in
%   maxima using 'latmax' and 'lmax'.
%
%   Compare to YLM which has a different (-1)^m phase factor.
%
%   Demo 2 is not working (the original demo in slepian_alpha is not
%   working either)
%   Demo 7 is not tested since I don't have the data (so if the original
%   demo is not working, this one won't work either)
%
% See also
%   XYZ2PLM, PLM2SPEC, TH2PL, PL2TH, YLM
%
% Special thanks to kwlewis@princeton.edu for spotting a bug.
% Last modified by
%   2024/08/13, williameclee@arizona.edu (@williameclee)
%   2023/11/20, fjsimons@alum.mit.edu (@fjsimons)

function varargout = plm2xyz(varargin)
    %% Initialisation
    % Demos
    if ischar(varargin{1}) || isstring(varargin{1})
        varargout = rundemos(varargin{:});

        if nargout == 0
            clear varargout
        end

        return
    end

    % Parse inputs
    [lmcosi, meshSize, c11cmn, lmax, latmax, Plm, beQuiet] = ...
        parseinputs(varargin{:});

    %% Computing the mesh
    % Lowest degree of the expansion
    lmin = lmcosi(1);
    % Highest degree (bandwidth of the expansion)
    L = lmcosi(end, 1);
    % Default resolution is the Nyquist degree; return equal sampling in
    % longitude and latitude; sqrt(L*(L+1)) is equivalent wavelength
    maxMeshSize = 180 / sqrt(L * (L + 1));
    meshSize = conddefval(meshSize, maxMeshSize);

    % When do you get a task bar?
    taskmax = 100;
    taskinf = 0;

    % But is it a grid or are they merely scattered points?
    if length(meshSize) == length(c11cmn)
        % It's a bunch of points!
        nLat = length(meshSize);
        nLon = length(c11cmn);
        % Colatitude vector in radians
        theta = deg2rad((90 - meshSize(:)'));
        % Longitude vector in radians
        phi = deg2rad(c11cmn(:)');

        % Initialize output vector
        r = zeros(nLat, 1);

        % Now if this is too large reduce lmax, our only recourse to hardcode
        ntb = 256;

        if round(sqrt(nLat)) >= ntb || round(sqrt(nLon)) >= ntb
            lmax = round(ntb);
        end

    elseif isscalar(meshSize) && length(c11cmn) == 4
        % It's a grid
        if meshSize > maxMeshSize && ~beQuiet
            disp('PLM2XYZ: You can do better! Ask for more spatial resolution')
            fprintf('Spatial sampling ALLOWED: %8.3f ; REQUESTED: %6.3f\n', ...
                maxMeshSize, meshSize)
        end

        % The number of longitude and latitude grid points that will be computed
        nLon = min(ceil((c11cmn(3) - c11cmn(1)) / meshSize + 1), latmax);
        nLat = min(ceil((c11cmn(2) - c11cmn(4)) / meshSize + 1), 2 * latmax + 1);

        % Initialize output grid
        r = zeros(nLat, nLon);

        % Longitude grid vector in radians
        phi = linspace(deg2rad(c11cmn(1)), deg2rad(c11cmn(3)), nLon);
        % Colatitude grid vector in radians
        theta = linspace(deg2rad(90 - c11cmn(2)), deg2rad(90 - c11cmn(4)), nLat);

        %    disp(sprintf('Creating %i by %i grid with resolution %8.3f',nlat,nlon,degres))
    else
        error('Make up your mind - is it a grid or a list of points?')
    end

    % Here we were going to build an option for a polar grid
    % But abandon this for now
    % [thetap,phip]=rottp(theta,phi,pi/2,pi/2,0);

    % Piecemeal degree ranges
    % Divide the degree range increments spaced such that the additional
    % number of degrees does not exceed addmup(lmax)
    % If this takes a long time, abort it
    els = zeros(1, 1);
    ind = 0;

    while els < L
        ind = ind + 1;
        % Take positive root
        els(ind + 1) = min(floor(max(roots( ...
            [1 3 -els(ind) ^ 2 - 3 * els(ind) - 2 * addmup(lmax)]))), L);

        if any(diff(els) == 0)
            error('Increase lmax as you are not making progress')
        end

    end

    % Now els contains the breakpoints of the degrees
    if ~all(diff(addmup(els)) <= addmup(lmax))
        error('The subdivision of the degree scale went awry')
    end

    % Here's the lspacings
    if length(els) > 2
        els = pauli(els, 2) + ...
            [0 0; ones(length(els) - 2, 1) zeros(length(els) - 2, 1)];
    end

    for ldeg = 1:size(els, 1)
        ldown = els(ldeg, 1);
        lup = els(ldeg, 2);
        % Construct the filename
        putputFolder = fullfile(getenv('IFILES'), 'LEGENDRE');
        outputFile = sprintf('LSSM-%i-%i-%i.mat', ldown, lup, nLat);
        outputPath = fullfile(putputFolder, outputFile);
        % ONLY COMPLETE LINEARLY SPACED SAMPLED VECTORS ARE TO BE SAVED!
        if (isfile(outputPath) && isequal(c11cmn, [0, 90, 360, -90])) ...
                && ~(size(els, 1) == 1 && ~isempty(Plm))
            % Get Legendre function values at linearly spaced intervals
            load(outputPath) %#ok<LOAD>

            if ~beQuiet
                fprintf('%s loaded %s\n', upper(mfilename), outputPath)
            end

            % AND TYPICALLY ANYTHING ELSE WOULD BE PRECOMPUTED, BUT THE GLOBAL
            % ONES CAN TOO! The Matlabpool check doesn't seem to work inside
        elseif size(els, 1) == 1 && ~isempty(Plm) && matlabpool('size') == 0 %#ok<*DPOOL>
            % disp(sprintf('Using precomputed workspace Legendre functions'))
        else
            % Evaluate Legendre polynomials at selected points
            try
                Plm = nan(length(theta), addmup(lup) - addmup(ldown - 1));
            catch ME
                error('\n %s \n\n Decrease lmax in PLM2XYZ \n', ME.message)
            end

            if lup - ldown > taskmax && length(meshSize) > taskinf
                h = waitbar(0, sprintf( ...
                    'Evaluating Legendre polynomials between %i and %i', ...
                    ldown, lup));
            end

            in1 = 0;
            in2 = ldown + 1;
            % Always start from the beginning in this array, regardless of lmin
            for l = ldown:lup

                % Never use Libbrecht algorithm... found out it wasn't that good
                Plm(:, in1 + 1:in2) = (legendre(l, cos(theta(:)'), 'sch') * sqrt(2 * l + 1))';

                in1 = in2;
                in2 = in1 + l + 2;

                if lup - ldown > taskmax && length(meshSize) > taskinf
                    waitbar((l - ldown + 1) / (lup - ldown + 1), h)
                end

            end

            if lup - ldown > taskmax && length(meshSize) > taskinf
                delete(h)
            end

            if length(c11cmn) == 4 && isequal(c11cmn, [0, 90, 360, -90])

                try
                    save(outputPath, 'Plm', '-v7.3')
                catch
                    save(outputPath, 'Plm')
                end

                if ~beQuiet
                    fprintf('%s saved %s\n', upper(mfilename), outputPath)
                end

            end

        end

        % Loop over the degrees
        more off
        %    disp(sprintf('PLM2XYZ Expansion from %i to %i',max(lmin,ldown),lup))
        for l = max(lmin, ldown):lup
            % Compute Schmidt-normalized Legendre functions at
            % the cosine of the colatitude (=sin(lat)) and
            % renormalize them to the area of the unit sphere

            % Remember the Plm vector always starts from ldown
            b = addmup(l - 1) + 1 - addmup(ldown - 1);
            e = addmup(l) - addmup(ldown - 1);

            plm = Plm(:, b:e)';

            m = 0:l;
            mphi = m(:) * phi(:)';

            % Normalization of the harmonics is to 4\ pi, the area of the unit
            % sphere: $\int_{0}^{\pi}\int_{0}^{2\pi}
            % |P_l^m(\cos\theta)\cos(m\phi)|^2\sin\theta\,d\theta d\phi=4\pi$.
            % Note the |cos(m\phi)|^2 d\phi contributes exactly \pi for m~=0
            % and 2\pi for m==0 which explains the absence of sqrt(2) there;
            % that fully normalized Legendre polynomials integrate to 1/2/pi
            % that regular Legendre polynomials integrate to 2/(2l+1),
            % Schmidt polynomials to 4/(2l+1) for m>0 and 2/(2l+1) for m==0,
            % and Schmidt*sqrt(2l+1) to 4 or 2. Note that for the integration of
            % the harmonics you get either 2pi for m==0 or pi+pi for cosine and
            % sine squared (cross terms drop out). This makes the fully
            % normalized spherical harmonics the only ones that consistently give
            % 1 for the normalization of the spherical harmonics.
            % Note this is using the cosines only; the "spherical harmonics" are
            % actually only semi-normalized.
            % Test normalization as follows (using inaccurate Simpson's rule):
            defval('tst', 0)

            if tst
                f = (plm .* plm)' .* repmat(sin(theta(:)), 1, l + 1);
                c = (cos(mphi) .* cos(mphi))';
                ntest = simpson(theta, f) .* simpson(phi, c) / 4 / pi;
                fprintf('Mean normalisation error l= %3.3i: %8.3e\n', ...
                    l, (sum(abs(1 - ntest))) / (l + 1))
                % For a decent test you would use "legendreprodint"
            end

            % Find the cosine and sine coefficients for this degree
            clm = shcos(lmcosi, l);
            slm = shsin(lmcosi, l);

            fac1 = repmat(clm, 1, nLon) .* cos(mphi) + ...
                repmat(slm, 1, nLon) .* sin(mphi);

            % Sum over all orders and (through loop) over all degrees
            if length(meshSize) == length(c11cmn)
                expa = sum(plm .* fac1, 1)';
                % Or diag(plm'*fac1) if you will
            elseif isscalar(meshSize) & length(c11cmn) == 4
                expa = plm' * fac1;
            end

            r = r + expa;
        end

    end

    lon = rad2deg(phi);
    lat = 90 - rad2deg(theta);

    % Prepare output
    varargout = {r, lon, lat, Plm, meshSize};
end

%% Subfunctions
function Outputs = rundemos(varargin)
    demoId = varargin{1};

    Outputs = cell(1);

    switch demoId
        case 'demo1'
            plm2xyz_demo1;
        case 'demo2'

            try
                fraction = varargin{2};
                plm2xyz_demo2(fraction);
            catch
                plm2xyz_demo2;
            end

        case 'demo3'
            plm2xyz_demo3;
        case 'demo4'
            plm2xyz_demo4;
        case 'demo5'
            plm2xyz_demo5;
        case 'demo6'
            Outputs = plm2xyz_demo6(nargout);
        case 'demo7'
            plm2xyz_demo7;
        case 'demo8'

            try
                meshSize = varargin{2};
                Outputs = plm2xyz_demo8(meshSize);
            catch
                Outputs = plm2xyz_demo8;
            end

    end

end

function varargout = parseinputs(varargin)
    meshSizeD = 1;
    c11cmnD = [0, 90, 360, -90];
    LmaxD = 720;
    latmaxD = Inf;
    p = inputParser;
    addRequired(p, 'lmcosi', @(x) isnumeric(x));
    addOptional(p, 'MeshSize', meshSizeD, ...
        @(x) (isnumeric(x) && isscalar(x)) || isempty(x));
    addOptional(p, 'c11cmn', c11cmnD, ...
        @(x) (isnumeric(x) && isvector(x)) || isempty(x));
    addOptional(p, 'Lmax', LmaxD, ...
        @(x) (isnumeric(x) && isscalar(x)) || isempty(x));
    addOptional(p, 'latmax', latmaxD, ...
        @(x) (isnumeric(x) && isscalar(x)) || isempty(x));
    addOptional(p, 'Plm', [], ...
        @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    lmcosi = p.Results.lmcosi;
    meshSize = conddefval(p.Results.MeshSize, meshSizeD);
    c11cmn = conddefval(p.Results.c11cmn, c11cmnD);
    lmax = conddefval(p.Results.Lmax, LmaxD);
    latmax = conddefval(p.Results.latmax, latmaxD);
    Plm = p.Results.Plm;
    beQuiet = logical(p.Results.BeQuiet);

    varargout = {lmcosi, meshSize, c11cmn, lmax, latmax, Plm, beQuiet};
end
