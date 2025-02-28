%% XYZ2PLM
% [lmcosi,dw,L2err]=XYZ2PLM(mesh,L,method,lat,lon,cnd)
%
% Forward real spherical harmonic transform in the 4pi normalized basis.
%
% Converts a spatially gridded field into spherical harmonics.
% For complete and regular spatial samplings [0 360 -90 90].
% If regularly spaced and complete, do not specify lat,lon.
% If not regularly spaced, mesh, lat and lon are column vectors.
%
% INPUT:
%
% mesh         Real-valued function whose transform we seek:
%               [1] MxN matrix of values corresponding to a regular (mesh)grid
%               defined by lat,lon as described below, OR
%               [2] an MNx1 vector of values corrsponding to a set of
%               latitude and longitude values given by lat,lon as below
% L             Bandwidth (maximum angular degree) or passband (two degrees)
%               Bandpass only implemented for method 'irr'
%               The maximum degree of the expansion is Nyquist checked
% method        'im'         By inversion (fast, accurate, preferred),
%                            uses FFT on equally spaced longitudes, ok to
%                            specify latitudes only as long as nat>=(L+1),
%                            note: works with the orthogonality of the
%                            cosine/sine of the longitude instead of with
%                            the orthogonality of the Legendre polynomials.
%               'gl'         By Gauss-Legendre integration (fast, inaccurate)
%                            note: resampling to GL integration points,
%                            uses FFT on equally spaced longitudes
%               'simpson'    By Simpson integation (fast, inaccurate),
%                            note: requires equidistant latitude spacing,
%                            uses FFT on equally spaced longitudes
%               'irr'        By inversion (irregular samplings)
%               'fib'        By Riemann sum on a Fibonacci grid (not done yet)
% lat           Latitude range for the grid or set of function values:
%               [1] if unspecified, we assume [90 -90] and a regular grid
%               [2] 1x2 vector [maximumlatitude minimumlatitude] in degrees
%               [3] an MNx1 vector of values with the explicit latitudes
% lon           Longitude range for the grid or set of function values:
%               [1] if unspecified, we assume [0 360] and a regular grid
%               [2] 1x2 vector [maximumlatitude minimumlatitude] in degrees
%               [3] an MNx1 vector of values with the explicit longitudes
% cnd           Eigenvalue tolerance in the irregular case
%
% OUTPUT:
%
% lmcosi        Matrix listing l,m,cosine and sine coefficients
% dw            Eigenvalue spectrum in the irregular case
% L2err         L2 error in the irregular case
%
% Note that the MEAN of the input data deviates from C(1), as sampled
% fields lose the orthogonality. The inversion approaches should recover
% the exact value of C(1), the true mean of the data, not the sample
% mean.
%
% lmcosi=xyz2plm(ones(randi(100),randi(100))); lmcosi(1,3) is close to one
%
% See also PLM2XYZ, PLM2SPEC, PLOTPLM, etc.
%
% Last modified by fjsimons-at-alum.mit.edu, 11/20/2023

function [lmcosi, dw, L2err] = xyz2plm_new(varargin)

    [fthph, L, method, lat, lon, cnd, beQuiet] = parseinputs(varargin{:});
    dw = [];
    L2err = [];

    t0 = clock;

    % Figure out if it's lowpass or bandpass
    lp = isscalar(L);
    bp = length(L) == 2;
    maxL = max(L);

    if bp
        minL = min(L);
    else
        minL = 0;
    end

    as = 0;
    % If no grid is specified, assumes equal spacing and complete grid
    if isempty(lat) && isempty(lon)
        % Test if data is 2D, and periodic over longitude... but not for 'irr'
        fthph = reduntest(fthph, beQuiet);
        polestest(fthph, beQuiet);
        % Make a complete grid
        nLon = size(fthph, 2);
        nLat = size(fthph, 1);
        % Nyquist wavelength
        Lnyq = min([ceil((nLon - 1) / 2) nLat - 1]);
        % Colatitude and its increment
        theta = linspace(0, pi, nLat);
        as = 1; % Equally spaced
        % Calculate latitude/longitude sampling interval; no wrap-around left
        % dtheta = pi / (nLat - 1);
        dphi = 2 * pi / nLon;

        switch method
                % Even without lat/lon can still choose the full inversion method
                % without Fourier transformation
            case 'irr'
                [LON, LAT] = meshgrid(linspace(0, 2 * pi * (1 - 1 / nLon), nLon), ...
                    linspace(pi / 2, -pi / 2, nLat));
                lat = LAT(:); lon = LON(:); fthph = fthph(:);
                theta = pi / 2 - lat;
                clear LON LAT
        end

    elseif isempty(lon)
        % If only latitudes are specified; make equal spacing longitude grid
        % Latitudes can be unequally spaced for 'im', 'irr' and 'gl'.
        fthph = reduntest(fthph, beQuiet);
        theta = (90 - lat) * pi / 180;
        dtheta = (lat(1) - lat(2)) * pi / 180;
        nLat = length(lat);
        nLon = size(fthph, 2);
        dphi = 2 * pi / nLon;
        Lnyq = min([ceil((nLon - 1) / 2) ceil(pi / dtheta)]);
    elseif length(lon) == 2 && length(lat) == 2 && ...
            lat(1) == 90 && lat(2) == -90 && diff(lon) == 360
        % This is a rotated complete map and requires fixing
        error('needs fixing')
        fthph = maprotate(fthph, [lon(1) lat(1) lon(2) lat(2)]);
        [lmcosi, dw] = xyz2plm(fthph, L, [], [], [], []);
        return
        % Is it a square? Force that to irregular method...
    else
        % Irregularly sampled data
        fthph = fthph(:);
        theta = deg2rad((90 - lat));
        lat = deg2rad(lat(:));
        lon = deg2rad(lon(:));
        nLon = length(lon);
        nLat = length(lat);
        % Nyquist wavelength
        adi = [abs(diff(sort(lon))); abs(diff(sort(lat)))];
        Lnyq = ceil(pi / min(adi(~ ~adi)));
        method = 'irr';
    end

    % Decide on the Nyquist frequency
    defval('maxL', Lnyq);
    % Never use Libbrecht algorithm... found out it wasn't that good
    defval('libb', 0)
    %disp(sprintf('Lnyq= %i ; expansion out to degree L= %i',Lnyq,L))

    if maxL > Lnyq || nLat < (maxL + 1)
        warning('XYZ2PLM: Function undersampled. Aliasing will occur.')
    end

    % Make cosine and sine matrices
    [m, l, mz] = addmon(maxL);
    lmcosi = [l m zeros(length(l), 2)];

    if bp
        lmcosi = lmcosi(addmup(L(1) - 1) + 1:end, :);
        % The l gets reused in the loop below but the m needs adapting
        m = m(addmup(L(1) - 1) + 1:end, :);
    end

    % Define evaluation points
    switch method
        case 'gl'
            % Highest degree of integrand will always be 2*maxL
            [w, x] = gausslegendrecof(2 * maxL, [], [-1 1]);
            % Function interpolated at Gauss-Legendre latitudes; 2D no help
            fthph = interp1(theta, fthph, acos(x), 'spline');
        case {'irr', 'simpson', 'im'}
            % Where else to evaluate the Legendre polynomials
            x = cos(theta);
        otherwise
            error('Specify valid method')
    end

    if lp
        fnpl = sprintf('%s/LSSM-%i-%i.mat', ...
            fullfile(getenv('IFILES'), 'LEGENDRE'), L, length(x));
    else
        fnpl = sprintf('%s/LSSM-%i-%i-%i.mat', ...
            fullfile(getenv('IFILES'), 'LEGENDRE'), L(1), L(2), length(x));
    end

    if exist(fnpl, 'file') == 2 & as == 1
        load(fnpl)
    else
        % Evaluate Legendre polynomials at selected points
        Plm = NaN(length(x), size(lmcosi, 1));

        if maxL > 200
            h = waitbar(0, 'Evaluating all Legendre polynomials');
        end

        in1 = 0;
        in2 = minL + 1;
        % Loop over the degrees
        for l = minL:maxL

            if libb == 0
                Plm(:, in1 + 1:in2) = ...
                    (legendre(l, x(:)', 'sch') * sqrt(2 * l + 1))';
            else
                Plm(:, in1 + 1:in2) = ...
                    (libbrecht(l, x(:)', 'sch') * sqrt(2 * l + 1))';
            end

            in1 = in2;
            in2 = in1 + l + 2;

            if L > 200
                waitbar((l + 1) / (L + 1), h)
            end

        end

        if maxL > 200
            delete(h)
        end

        if as == 1
            save(fnpl, 'Plm')
        end

    end

    switch method
        case 'irr'
            Plm = [Plm .* cos(lon(:) * m(:)'), Plm .* sin(lon(:) * m(:)')];
            % Add these into the sensitivity matrix
            [C, ~, ~, ~, L2err, ~, dw] = datafit(Plm, fthph, [], [], cnd);
            lmcosi(:, 3) = C(1:end / 2);
            lmcosi(:, 4) = C(end / 2 + 1:end);
        case {'im', 'gl', 'simpson'}
            % Perhaps demean the data for Fourier transform
            defval('dem', 0)

            if dem
                meanm = mean(fthph, 2);
                fthph = fthph - repmat(meanm, 1, nLon);
            end

            % Calculate integration over phi by the fast Fourier
            % transform. Integration of real input field with respect to the second
            % dimension of r, at  wavenumber m, thus at constant latitude. You get
            % as many wavenumbers m as there are longitudes; only use to L. With
            % Matlab's FFT, need to multiply by sampling interval.
            gfft = dphi * fft(fthph, nLon, 2);

            if dem
                % Add the azimuthal mean back in there
                gfft(:, 1) = 2 * pi * meanm;
            end

            % Note these things are only half unique - the maximum m is nlon/2
            % But no Nyquist theory exists for the Legendre transform...
            a = real(gfft);
            b = -imag(gfft);
            in1 = 0;
            in2 = minL + 1;
        otherwise
            error('Specify valid method')
    end

    if bp
        error('needs fixing')
    end

    switch method
        case 'im'
            % Loop over the orders. This speeds it up versus 'irr'
            for ord = 0:L
                a(:, 1) = a(:, 1) / 2;
                b(:, 1) = b(:, 1) / 2;
                Pm = Plm(:, mz(ord + 1:end) + ord) * pi;
                [lmcosi(mz(ord + 1:end) + ord, 3)] = ...
                    datafit(Pm, a(:, ord + 1), [], [], cnd);
                [lmcosi(mz(ord + 1:end) + ord, 4)] = ...
                    datafit(Pm, b(:, ord + 1), [], [], cnd);
            end

        case 'simpson'
            % Loop over the degrees. Could go up to l=nlon if you want
            for l = 0:L
                % Integrate over theta using Simpson's rule
                clm = simpson(theta, ...
                    repmat(sin(theta(:)), 1, l + 1) .* a(:, 1:l + 1) .* Plm(:, in1 + 1:in2));
                slm = simpson(theta, ...
                    repmat(sin(theta(:)), 1, l + 1) .* b(:, 1:l + 1) .* Plm(:, in1 + 1: ...
                    in2));
                in1 = in2;
                in2 = in1 + l + 2;
                % And stick it in a matrix [l m Ccos Csin]
                lmcosi(addmup(l - 1) + 1:addmup(l), 3) = ...
                    clm(:) / 4 / pi;
                lmcosi(addmup(l - 1) + 1:addmup(l), 4) = ...
                    slm(:) / 4 / pi;
            end

        case 'gl'
            % Loop over the degrees. Could go up to l=nlon if you want
            for l = 0:L
                % Integrate over theta using Gauss-Legendre integration
                clm = sum(a(:, 1:l + 1) .* (diag(w) * Plm(:, in1 + 1:in2)));
                slm = sum(b(:, 1:l + 1) .* (diag(w) * Plm(:, in1 + 1:in2)));
                in1 = in2;
                in2 = in1 + l + 2;
                % And stick it in a matrix [l m Ccos Csin]
                lmcosi(addmup(l - 1) + 1:addmup(l), 3) = ...
                    clm(:) / 4 / pi;
                lmcosi(addmup(l - 1) + 1:addmup(l), 4) = ...
                    slm(:) / 4 / pi;
            end

            % rnk = [];
    end

    % Get rid of machine precision error
    lmcosi(abs(lmcosi(:, 3)) < eps, 3) = 0;
    lmcosi(abs(lmcosi(:, 4)) < eps, 4) = 0;

    % fprintf('XYZ2PLM (Analysis using %s)  took %8.4f s\n', method, etime(clock, t0))
end

%% Subfunctions
function varargout = parseinputs(varargin)
    LD = 18;
    methodD = 'im';
    p = inputParser;
    addRequired(p, 'fthph', ...
        @(x) isnumeric(x) && ismatrix(x));
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) && (isscalar(x) || length(x) == 2));
    addOptional(p, 'method', methodD, ...
        @(x) ischar(x) || isstring(x) && ismember(x, {'im', 'gl', 'simpson', 'irr'}));
    addOptional(p, 'lat', [], ...
        @(x) isnumeric(x));
    addOptional(p, 'lon', [], ...
        @(x) isnumeric(x));
    addOptional(p, 'cnd', [], ...
        @(x) true);
    addParameter(p, 'BeQuiet', false, ...
        @istruefalse);
    parse(p, varargin{:});

    fthph = p.Results.fthph;
    L = conddefval(p.Results.L, LD);
    method = conddefval(p.Results.method, methodD);
    lat = p.Results.lat;
    lon = p.Results.lon;
    cnd = p.Results.cnd;
    beQuiet = logical(p.Results.BeQuiet);

    varargout = {fthph, L, method, lat, lon, cnd, beQuiet};
end

function grd = reduntest(grd, beQuiet)
    % Tests if last longitude repeats last (0,360)
    % and removes last data column
    if sum(abs(grd(:, 1) - grd(:, end))) >= size(grd, 2) * eps * 10 && ~beQuiet
        fprintf('%s: Data violate wrap-around by %8.4e\n', ...
            upper(mfilename), sum(abs(grd(:, 1) - grd(:, end))))
    end

    grd = grd(:, 1:end - 1);
end

function polestest(grd, beQuiet)
    % Tests if poles (-90,90) are identical over longitudes
    var1 = var(grd(1, :));
    var2 = var(grd(end, :));

    if var1 > eps * 10 | var2 > eps * 10 && ~beQuiet
        fprintf('%s: Poles violated by %8.4e and %8.4e\n', upper(mfilename), var1, var2)
    end

end
