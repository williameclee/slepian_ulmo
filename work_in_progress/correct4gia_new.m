% [thedates,GIAt,GIAtU,GIAtL,trend]=CORRECT4GIA(thedates,model,TH,L)
%
% This function accepts an array of dates (in matlab datenum format) and
% calculates the surface mass density change from a certain GIA model at
% each date with reference to the first date.  That is, the GIA models are
% saved as yearly rates of surface mass density change, so there will be a
% different field for each date.
%
% If you want lmcosi fields for GIA, then give only 2 arguments.  If you
% want the GIA projected into a Slepian basis, give that basis.
%
%
% INPUT:
%
% thedates     An array of dates in matlab DATENUM format.  The first date
%                is used as the reference date. [default: monthly points
%                during 2004]
% model        Which GIA model you want.  See references. Options include:
%                'Paulson07'    A model based on the ICE-5G ice load model
%                               of Peltier (2004).  Suitable for both
%                               Antarctica and Greenland.  As corrected by
%                               Geruo A and J. Wahr.  Global.
%                'Wangetal08'   A model based on the older ICE-4G ice
%                               model, and viscosity which varies
%                               laterally.  Suitable for Greenland.
%                'IJ05_R2'      A model based on the Ivins et al. (2013)
%                               ice model.  Suitable for Antarctica.
%                'IJ05'         A model based on the Ivins and James (2005)
%                               ice model.  Suitable for Antarctica.
%                'W12a_v1'      A "best" model from Whitehouse et al (2012)
%                               Suitable only for Antarctica.
%
% TH         Optional [default nothing]:Radius of the concentration
%              region (degrees) OR
%              'england', 'eurasia',  'namerica', 'australia', 'greenland'
%              'africa', 'samerica', 'amazon', 'orinoco', in which case
%              you must have phi,theta,omega all equal to zero OR
%              [lon lat] an ordered list defining a closed curve [degrees]
%              OR a cell containing a region and a buffer such
%              as {'greenland' 0.5}
% Lwindow    Optional [default nothing]: Bandwidth of the window [default
%              if you give a region: bandwidth of the data]
%          If you gave a concentration radius for TH, then you need to
%          specify these:
% phi        Longitude of the center of the tapers (degrees)
% theta      Colatitude of the center of the tapers (degrees)
% omega      Anticlockwise azimuthal rotation of the tapers (degrees)
%
%
%
% OUTPUT:
%
% GIAt          If you want SH coefficients for the GIA, then this will be
%                a 3D matrix of GIA fields for each date requested.  The
%                first dimension are the dates.  The second and third
%                dimensions are the familiar lmcosi dimensions
%                (i.e. (L+1)^2 by 4)  Units are kg/m^2
%                OR if you asked for a Slepian projection of this data,
%                then it will be a matrix like "slept" where the first
%                dimension is time and the second dimension is
%                Slepian coefficient.
% thedates      Your dates back to you.
% GIAtU         Same as GIAt, but if the model has an upper bound
% GIAtL         Same as GIAt, but if the model has a lower bound
% trend         This is the magnitude of the correction in units of Gt/yr.
%                This output will only work if you gave a Slepian basis.
%
%
% NOTES:  It is left to the user to obtain the aforementioned GIA models
% from the authors and save them as appropriate mat files (this function
% assumes rate of change of surface mass density per year).  Perhaps in
% the future these models will be distributed collectively.
%
% SEE ALSO:  PLM2POT
% REFERENCES:
%   Paulson, A., S. Zhong, and J. Wahr. Inference of mantle viscosity from
%    GRACE and relative sea level data, Geophys. J. Int. (2007) 171,
%    497–508. doi: 10.1111/j.1365-246X.2007.03556.x
%
%   Geruo, A., Wahr, J. & Zhong, S. Computations of the viscoelastic response of a
%    3-D compressible Earth to surface loading: An application to Glacial Isostatic
%    Adjustment in Antarctica and Canada. Geophys. J. Int. 192, 557-572 (2013).
%
%   Ivins, E. R., T. S. James, J. Wahr, E. J. O. Schrama, F. W. Landerer,
%   and K. M. Simon. Antarctic contribution to sea level rise observed by
%   GRACE with improved GIA correction, Journal of Geophysical Research:
%   Solid Earth, vol. 118, 3126-3141, doi: 10.1002/jgrb.50208, 2013.
%
%   Ivins, E. R., T. S. James. Antarctic glacial isostatic adjustment: a
%   new assessment, Antarctic Science, vol. 17(4), 541-553,
%   doi: 10.1017/S0954102005002968, 2005.
%
%   Wang, H., P. Wu, and W. van der Wal. Using postglacial sea level, crustal
%    velocities and gravity-rate-of-change to constrain the influence of
%    thermal effects on mantle lateral heterogeneities, Journal of
%    Geodynamics (2008) 46, 104-117. doi: 10.1016/j.jog.2008.03.003
%
%   Whitehouse, P.L., Bentley, M.J., Milne, G.A., King, M.A.,
%    Thomas, I.D., 2012. A new glacial isostatic adjustment model for
%    Antarctica: calibrated and tested using observations of relative
%    sea-level change and present-day uplift rates. Geophysical Journal
%    International 190, 1464-1482. doi:10.1111/j.1365-246X.2012.05557.x
%
% Last modified by
%   charig-at-email.arizona.edu on 5/23/2017
%   williameclee-at-arizona.edu on 7/12/2024

function varargout = correct4gia_new(varargin)
    %% Initialisation
    % Parse inputs
    [time, model, domain, L, phi, theta, omega, ~, beQuiet] = parseinputs(varargin{:});

    % Where the model save files are kept
    if strncmp(model, 'Morrow', 6)
        inputFolder = fullfile(getenv('IFILES'), 'GIA', model(1:6));
    else
        inputFolder = fullfile(getenv('IFILES'), 'GIA', model);
    end

    % And the appropriate name
    inputPath = fullfile(inputFolder, sprintf('%s_SD.mat', model));

    % Load this data (saved as lmcosiM)
    load(inputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');

    %% Main
    % Convert the model from change per year to change per day
    lmcosiM(:, 3:4) = lmcosiM(:, 3:4) / 365;

    if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
        lmcosiU(:, 3:4) = lmcosiU(:, 3:4) / 365;
        lmcosiL(:, 3:4) = lmcosiL(:, 3:4) / 365;
    end

    % Reference the date string to the first date
    if ~isscalar(time)
        deltatime = time - time(1);
    else
        deltatime = 365;
    end

    % If we have a Slepian basis, do that, if not just do plms
    if exist('domain', 'var')

        % Figure out if it's lowpass or bandpass
        % lp = length(L) == 1;
        % bp = length(L) == 2;
        maxL = max(L);
        % The spherical harmonic dimension
        % ldim = (L(2 - lp) + 1) ^ 2 - bp * L(1) ^ 2;

        % Project the model
        [~, ~, ~, lmcosiW, ~, ~, ~, ~, ~, ronm] = addmon(maxL);

        if isa(domain, 'GeoDomain') || ismatrix(domain)
            [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, "BeQuiet", beQuiet);
        else
            [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, phi, theta, omega, "BeQuiet", beQuiet);
        end

        % if isnumeric(domain)
        %     [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, phi, theta, omega, "BeQuiet", beQuiet);
        % elseif iscell(domain) && length(domain) >= 3
        %     [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, "BeQuiet", beQuiet);
        % else
        %     [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, "MoreRegionSpecs", moreDomainSpecs, "BeQuiet", beQuiet);
        % end

        if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
            falphaU = plm2slep_new(lmcosiU, domain, L, "BeQuiet", beQuiet);
            falphaL = plm2slep_new(lmcosiL, domain, L, "BeQuiet", beQuiet);
        end

        % Scale to the new dates
        GIAt = zeros([length(deltatime), length(falpha)]);

        if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
            GIAtU = zeros([length(deltatime), length(falphaU)]);
            GIAtL = zeros([length(deltatime), length(falphaL)]);
        end

        for i = 1:length(deltatime)
            GIAt(i, :) = falpha * deltatime(i);

            if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
                GIAtU(i, :) = falphaU * deltatime(i);
                GIAtL(i, :) = falphaL * deltatime(i);
            end

        end

        % How large is this signal we just made?
        % To answer this we integrate the basis functions and multiply by the
        % coefficients representing the annual rate.
        % Note: the falpha from the models should already be in units of
        % surface mass density change per year.
        for j = 1:round(N)
            cosi = lmcosiW(:, 3:4);
            cosi(ronm) = G(:, j);
            CC{j} = [lmcosiW(:, 1:2), cosi];
        end

        if isa(domain, 'GeoDomain') || ismatrix(domain)
            eigfunINT = integratebasis_new(CC, domain, round(N));
        else
            eigfunINT = integratebasis_new(CC, domain, round(N), phi, theta);
        end

        % Since Int should have units of (fn * m^2), need to go from fractional
        % sphere area to real area.  If the fn is surface density, this output is
        % in kilograms.  Then change the units from kg to Gt in METRIC tons
        eigfunINT = eigfunINT * (4 * pi * 6370e3 ^ 2) / 1e12;
        % functionintegrals = eigfunINT;

        % Now multiply by the appropriate slepcoffs to get the months
        % This becomes alpha by months
        %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
        %functimeseries = sleptdelta(:,1:N)';

        total = eigfunINT(1:round(N)) .* (falpha(1:round(N)) * 365); % Back to per year

    else % Just do the plms
        GIAt = zeros([length(deltatime), size(lmcosiM)]);

        if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
            GIAtU = zeros([length(deltatime), size(lmcosiU)]);
            GIAtL = zeros([length(deltatime), size(lmcosiL)]);
        end

        for i = 1:length(deltatime)
            GIAt(i, :, 3:4) = lmcosiM(:, 3:4) * deltatime(i);

            if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
                GIAtU(i, :, 1:2) = lmcosiU(:, 1:2);
                GIAtU(i, :, 3:4) = lmcosiU(:, 3:4) * deltatime(i);
                GIAtL(i, :, 1:2) = lmcosiL(:, 1:2);
                GIAtL(i, :, 3:4) = lmcosiL(:, 3:4) * deltatime(i);
            end

        end

        total = 0;
    end

    if size(GIAt, 1) == 1
        GIAt = squeeze(GIAt);
    end

    %% Returning requested outputs
    if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
        varargout = {GIAt, time, GIAtU, GIAtL};
    else
        varargout = {GIAt, time, [], [], total};
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    timeD = datenum(2004, 1:12, 1); %#ok<DATNM>
    modelD = 'Paulson07';
    domainD = {'greenland', 0.5};
    LD = 60;
    phiD = 0;
    thetaD = 0;
    omegaD = 0;
    p = inputParser;
    addOptional(p, 'Time', timeD, ...
        @(x) isnumeric(x) || isdatetime(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || isempty(x));
    addOptional(p, 'Domain', domainD, ...
        @(x) isnumeric(x) || iscell(x) || ischar(x) || isa(x, 'GeoDomain') || isempty(x));
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'phi', phiD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'theta', thetaD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'omega', omegaD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'MoreDomainSpecs', {}, @iscell);
    addParameter(p, 'BeQuiet', false, @islogical);

    parse(p, varargin{:});
    time = conddefval(p.Results.Time, timeD);
    model = conddefval(p.Results.Model, modelD);
    domain = conddefval(p.Results.Domain, domainD);
    L = conddefval(p.Results.L, LD);
    phi = conddefval(p.Results.phi, phiD);
    theta = conddefval(p.Results.theta, thetaD);
    omega = conddefval(p.Results.omega, omegaD);
    moreDomainSpecs = p.Results.MoreDomainSpecs;
    beQuiet = p.Results.BeQuiet;

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, moreDomainSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2}, moreDomainSpecs{:});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:}, moreDomainSpecs{:});
    end

    if isdatetime(time)
        time = datenum(time); %#ok<DATNM>
    end

    varargout = {time, model, domain, L, phi, theta, omega, moreDomainSpecs, beQuiet};

end
