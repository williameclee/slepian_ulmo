% [dates, GIAt, GIAtU, GIAtL, trend] = correct4gia(dates, model, domain, L)
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
% Input arguments
%   thedates - An array of dates in matlab DATENUM format.  The first date
%                is used as the reference date. [default: monthly points
%                during 2004]
%   model - Which GIA model you want.  See references. Options include:
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
%   domain - Optional [default nothing]:Radius of the concentration
%              region (degrees) OR
%              'england', 'eurasia',  'namerica', 'australia', 'greenland'
%              'africa', 'samerica', 'amazon', 'orinoco', in which case
%              you must have phi,theta,omega all equal to zero OR
%              [lon lat] an ordered list defining a closed curve [degrees]
%              OR a cell containing a region and a buffer such
%              as {'greenland' 0.5}
%   L - Optional [default nothing]: Bandwidth of the window [default
%              if you give a region: bandwidth of the data]
%          If you gave a concentration radius for TH, then you need to
%          specify these:
%   phi - Longitude of the center of the tapers (degrees)
%   theta - Colatitude of the center of the tapers (degrees)
%   omega - Anticlockwise azimuthal rotation of the tapers (degrees)
%
% Output arguments
%   GIAt - If you want SH coefficients for the GIA, then this will be
%                a 3D matrix of GIA fields for each date requested.  The
%                first dimension are the dates.  The second and third
%                dimensions are the familiar lmcosi dimensions
%                (i.e. (L+1)^2 by 4)  Units are kg/m^2
%                OR if you asked for a Slepian projection of this data,
%                then it will be a matrix like "slept" where the first
%                dimension is time and the second dimension is
%                Slepian coefficient.
%   dates - Your dates back to you.
%   GIAtU - Same as GIAt, but if the model has an upper bound
%   GIAtL - Same as GIAt, but if the model has a lower bound
%   trend - This is the magnitude of the correction in units of Gt/yr.
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
    [time, model, domain, L, phi, theta, omega, ~, beQuiet] = ...
        parseinputs(varargin{:});

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
    hasBounds = exist('lmcosiU', 'var') && exist('lmcosiL', 'var');

    %% Main
    % Reference the date string to the first date
    if ~isscalar(time)
        deltaYear = (time - time(1)) / 365;
    else
        deltaYear = 1;
    end

    % If we have a Slepian basis, do that, if not just do plms
    if exist('domain', 'var')
        % Project the model
        if isa(domain, 'GeoDomain') || ismatrix(domain)
            [slep, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, ...
                "BeQuiet", beQuiet);
        else
            [slep, ~, N, ~, G] = plm2slep_new(lmcosiM, domain, L, ...
                phi, theta, omega, "BeQuiet", beQuiet);
        end

        if hasBounds
            slepU = plm2slep_new(lmcosiU, domain, L, "BeQuiet", beQuiet);
            slepL = plm2slep_new(lmcosiL, domain, L, "BeQuiet", beQuiet);
        end

        % Scale to the new dates
        slept = zeros([length(deltaYear), length(slep)]);

        if hasBounds
            sleptU = zeros([length(deltaYear), length(slepU)]);
            sleptL = zeros([length(deltaYear), length(slepL)]);
        end

        for iTime = 1:length(deltaYear)
            slept(iTime, :) = slep * deltaYear(iTime);

            if hasBounds
                sleptU(iTime, :) = slepU * deltaYear(iTime);
                sleptL(iTime, :) = slepL * deltaYear(iTime);
            end

        end

        % How large is this signal we just made?
        % To answer this we integrate the basis functions and multiply by the
        % coefficients representing the annual rate.
        % Note: the falpha from the models should already be in units of
        % surface mass density change per year.
        J = round(N);

        if isa(domain, 'GeoDomain') || ismatrix(domain)
            eigfunINT = integratebasis_new(G, domain, J, ...
                "BeQuiet", beQuiet);
        else
            eigfunINT = integratebasis_new(G, domain, J, phi, theta, ...
                "BeQuiet", beQuiet);
        end

        % Since Int should have units of (fn * m^2), need to go from fractional
        % sphere area to real area.  If the fn is surface density, this output is
        % in kilograms.  Then change the units from kg to Gt in METRIC tons
        eigfunINT = eigfunINT * (4 * pi * 6370e3 ^ 2);

        % Now multiply by the appropriate slepcoffs to get the months
        % This becomes alpha by months
        total = slept(:, 1:J) * eigfunINT(1:J)';
        total = total / 1e12; % Convert to gigaton

    else % Just do the plms
        slept = zeros([length(deltaYear), size(lmcosiM)]);

        if hasBounds
            sleptU = zeros([length(deltaYear), size(lmcosiU)]);
            sleptL = zeros([length(deltaYear), size(lmcosiL)]);
        end

        for iTime = 1:length(deltaYear)
            slept(iTime, :, 3:4) = lmcosiM(:, 3:4) * deltaYear(iTime);

            if hasBounds
                sleptU(iTime, :, 1:2) = lmcosiU(:, 1:2);
                sleptU(iTime, :, 3:4) = lmcosiU(:, 3:4) * deltaYear(iTime);
                sleptL(iTime, :, 1:2) = lmcosiL(:, 1:2);
                sleptL(iTime, :, 3:4) = lmcosiL(:, 3:4) * deltaYear(iTime);
            end

        end

        total = 0;
    end

    %% Returning requested outputs
    slept = squeeze(slept);
    if hasBounds
        sleptU = squeeze(sleptU);
        sleptL = squeeze(sleptL);
        varargout = {slept, time, sleptU, sleptL};
    else
        varargout = {slept, time, [], [], total};
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
        @(x) isnumeric(x) || iscell(x) || ischar(x) || ...
        isa(x, 'GeoDomain') || isempty(x));
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

    varargout = {time, model, domain, L, phi, theta, omega, ...
        moreDomainSpecs, beQuiet};

end
