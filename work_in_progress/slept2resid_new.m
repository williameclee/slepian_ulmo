% [sleptCoeffsSignal, sleptCoeffsResid, ftests, extravalues, total, ...
%   alphavarall, totalparams, totalparamerrors, totalfit, ...
%   functionintegrals, alphavar]
%   = SLEPT2RESID(sleptCoeffs, dates, fitwhat, givenerrors, ...
%       specialterms, CC, TH)
%
% Takes a time series of Slepian coefficients and fits a desired
% combination of functions (e.g. secular, annual, semiannual, etc.) to
% each coefficient in order to separate "signal" from residual.
%
% You can choose to fit either a mean, linear, quadratic, or cubic fuction
% as your "base" function to each Slepian coefficient, by using "fitwhat".
% In these cases of higher functions, they are only used if they pass an
% F-test for variance reduction.  For example, a cubic term is only used if
% it is significantly better than a quadratic function.
%
% If you also provide the Slepian functions (CC) and region (TH) for this
% localization, then we assume that you want the fitting to the integrated
% functions. i.e. If you have surface density, then using the integrated
% Slepian function would give you total mass change in a region.  In
% addition, there will be a fitting of the combination (total) of the
% Slepian functions up to the Shannon number for this localization.
%
% INPUT:
%
% sleptCoeffs   The time series of Slepian coefficients. This should be a
%               two dimensional matrix (not a cell array), where the first
%               dimension is time, and the second dimension are Slepian
%               coefficients sorted by global eigenvalue.
% dates         An array of dates corresponding to the sleptCoeffs timeseries.
%               The unit is datenum.
%               ***It is assumed that 'dates' matches 'sleptCoeffs'. If
%               'dates' is longer than 'sleptCoeffs' then it is assumed that
%               these are extra dates that you would like sleptCoeffsSignal
%               to be evaluated at. This is useful if you want to
%               interpolate to one of the GRACE missing months that is
%               important, such as Jan 2011.
% fitwhat       The functions that you would like to fit to the time series
%               data. The format of this array is as follows:
%               [order periodic1 periodic2 periodic3 etc...] where
%               - order is either 0/1/2/3 to fit UP TO either a
%                   mean/linear/quadratic/cubic function (if increasing
%                   the power of the function reduces variance enough)
%               - periodic1 is the period in days of a function (i.e.
%                   365.0)
%               Any # of desired periodic functions [days] (including zero)
%               can be included.
% givenerrors   These are given errors, if you have them. In this case a
%               weighted inversion is performed. givenerrors should be the
%               same dimensions of sleptCoeffs.
% specialterms  A cell array such as {2 'periodic' 1460}. At the moment
%               this is pretty specific to our needs, but could be
%               expanded later.
% CC            A cell array of the localization Slepian functions
% TH            The region (proper string or XY coordinates) that you did
%               the localization on (so we can integrate)
% N             Number of largest eigenfunctions in which to expand. By
%               default rounds to the Shannon number.
%
% OUTPUT:
%
% sleptCoeffsSignal The least-squares fitted function for each Slepian
%               coefficient evaluated at those months in the same format
% sleptCoeffsResid  Residual time series for each Slepian coefficients,
%               ordered as they were given, presumably by eigenvalue
%               [nmonths x (Lwindow+1)^2]
% ftests        A matrix, such as [0 1 1] for each Slepian coefficient,
%               on whether the fits you requested passed an F-test for
%               significance.
% extravalues   These are the values of sleptCoeffsSignal evaluated at your
%               extra dates which were tacked onto 'dates'
% total         The time series of the combined mass change from N Slepian
%               functions (i.e. the combined data points)
%               The unit is Gt.
% alphavarall   The time averaged variance on each data point (from error
%               propogation from each individual function). These values
%               are the same for every point.
% totalparams   The parameters of the fits to the total. This is a 4-by-j
%               matrix where j are the fits you wanted. Zeros fill out the
%               unused parameters. For example, if you want just a 2nd
%               order polynomial, you will get [intercept intercept; slope
%               slope; 0 quadratic; 0 0]
%               The unit is Gt/day^p.
% totalparamerrors  The 95% confidence value on this slope. This is
%               computed using the critical value for a t-distribution.
%               At the moment, totalparams and totalparamerrors and just
%               the values for a linear fit. Sometime later maybe change
%               this to be potentially a quadratic fit as well.
%               The unit is Gt/yr.
% totalfit      Datapoints for the best fit line, so you can plot it, and
%               datapoints for the confidence intervals of this linear fit,
%               so you can plot them.  NOTE: To plot these you should use
%               totalfit(:,2)+totalfit(:,3) and totalfit(:,2)-totalfit(:,3)
% functionintegrals This is a vector of the integrals of the Slepian
%               functions, up to the Shannon number. With this we can
%               multiply by the data or sleptCoeffsSignal to get the
%               changes in total mass over time for each function.
%               The unit is Gt.
% alphavar      The variances on the individual alpha function time series
%               (INTEGRALS). This is calculated by making a covariance
%               matrix from sleptCoeffsResid and then doing the matrix
%               multiplication with the eigenfunction integrals.
%
% Last modified by
%   charig-at-princeton.edu  6/26/2012
%   williameclee-at-arizona.edu  6/07/2024

function varargout = slept2resid_new(varargin)
    % Add path to the auxiliary function(s)
    addpath(fullfile(fileparts(mfilename('fullpath')), 'aux'));
    %% Initialisation and preallocation
    sleptCoeffsDefault = 'grace2slept(''CSR'',''greenland'',0.5,60,[],[],[],[],''SD'')';
    datesDefault = datenum(2004, 1:12, 1); %#ok<DATNM>
    fitwhatDefault = [3, 365.0];
    givenerrorsDefault = [];
    % givenerrorsDefault = ones(size(sleptCoeffs));
    specialtermsDefault = {NaN};
    CCDefault = [];
    THDefault = [];
    NDefault = [];

    p = inputParser;
    addOptional(p, 'sleptCoeffs', sleptCoeffsDefault);
    addOptional(p, 'dates', datesDefault, @(x) isnumeric(x) || isdatetime(x));
    addOptional(p, 'Fit', fitwhatDefault);
    addOptional(p, 'givenerrors', givenerrorsDefault);
    addOptional(p, 'specialterms', specialtermsDefault);
    addOptional(p, 'CC', CCDefault);
    addOptional(p, 'TH', THDefault, ...
        @(x) ischar(x) || isnumeric(x) || iscell(x) || isempty(x));
    addOptional(p, 'N', NDefault);
    addOptional(p, 'MoreRegionSpecs', {});
    addParameter(p, 'OutputUnit', 'original', @(x) ischar(x) && ismember(x, {'original', 'year'}));
    parse(p, varargin{:});
    sleptCoeffs = conddefval(p.Results.sleptCoeffs, sleptCoeffsDefault);
    dates = conddefval(p.Results.dates, datesDefault);
    fitwhat = conddefval(p.Results.Fit, fitwhatDefault);
    givenerrors = conddefval(p.Results.givenerrors, givenerrorsDefault);
    specialterms = conddefval(p.Results.specialterms, specialtermsDefault);
    CC = conddefval(p.Results.CC, CCDefault);
    TH = conddefval(p.Results.TH, THDefault);
    N = conddefval(p.Results.N, NDefault);
    moreRegionSpecs = p.Results.MoreRegionSpecs;
    OutputUnit = p.Results.OutputUnit;

    if isdatetime(dates)
        dates = datenum(dates); %#ok<DATNM>
    end

    if iscell(TH) && length(TH) > 2
        moreRegionSpecs = {TH{3:end}, moreRegionSpecs{:}}; %#ok<CCAT>
    end

    if ischar(sleptCoeffs)
        % Evaluate the specified expression
        [sleptCoeffs, ~, dates, TH, ~, CC] = eval(sleptCoeffs);
    else
        givenerrors = conddefval(givenerrors, ones(size(sleptCoeffs)));
    end

    % Dates
    dates = dates(:)';
    nMonth = length(dates);

    if nMonth == size(sleptCoeffs, 1) % do nothing
        dateExtra = [];
        moredates = false;
    elseif nMonth > size(sleptCoeffs, 1) % pull off extra dates
        dateExtra = dates((size(sleptCoeffs, 1) + 1):end);
        dates = dates(1:size(sleptCoeffs, 1));
        moredates = true;
    else
        error('Is dates shorter than slept?')
    end

    % We will do a scaling to improve the solution
    dateMean = mean(dates);
    dateStd = std(dates);
    dateNormalised = (dates - dateMean) / dateStd;
    dateExtraNormalised = (dateExtra - dateMean) / dateStd;

    % Periodic components
    nOmega = length(fitwhat(2:end));

    if nOmega > 0
        % The (not angular) frequencies being fitted in [1/days]
        freq = 1 ./ fitwhat(2:end);
        % Rescale these to our new xprime
        freq = freq * dateStd;
    else
        freq = [];
    end

    % Slepian coefficients
    nSleptCoeff = size(sleptCoeffs, 2);
    % Preallocate the residuals and the evaluated fitted function set
    sleptCoeffsSignal = zeros(size(sleptCoeffs));
    sleptCoeffsResid = zeros(size(sleptCoeffs));
    extravalues = zeros([length(dateExtra), nSleptCoeff]);
    ftests = zeros([nSleptCoeff, fitwhat(1)]);
    % Figure out the orders and degrees of this setup
    % BUT orders and degrees have lost their meaning since slept should be
    % ordered by eigenvalue

    %% G matrix assembly
    % We will have the same number of G matrices as order of polynomial
    % fit. These matrices are smallish, so make all 3 regardless of whether
    % you want them all.
    G1 = []; % For line fits
    G2 = []; % For quadratic fits
    G3 = []; % For cubic fits
    % Mean term
    if fitwhat(1) >= 0
        G1 = [G1, ones(size(dateNormalised'))];
        G2 = [G2, ones(size(dateNormalised'))];
        G3 = [G3, ones(size(dateNormalised'))];
    end

    % Secular term
    if fitwhat(1) >= 1
        G1 = [G1, dateNormalised'];
        G2 = [G2, dateNormalised'];
        G3 = [G3, dateNormalised'];
    end

    % Quadratic term
    if fitwhat(1) >= 2
        G2 = [G2, dateNormalised' .^ 2];
        G3 = [G3, dateNormalised' .^ 2];
    end

    % Cubic term
    if fitwhat(1) == 3
        G3 = [G3, dateNormalised' .^ 3];
    end

    % Periodic terms
    if nOmega > 0
        % Angular frequency in radians/(rescaled day) of the periodic terms
        phase = repmat(freq, nMonth, 1) * 2 * pi ... % change into angular
            .* repmat(dateNormalised', 1, nOmega);
        G1 = [G1, cos(phase), sin(phase)];
        G2 = [G2, cos(phase), sin(phase)];
        G3 = [G3, cos(phase), sin(phase)];

        % Create our specialterms G if we have it. At the moment this is
        % just an additional periodic function, but in the future we could
        % add something else.
        if ~isnan(specialterms{1})
            % Strip off the previous periodic terms here and REPLACE with
            % freqSpec
            freqSpec = [freq, dateStd / specialterms{3}];
            phaseSpec = repmat(freqSpec, nMonth, 1) * 2 * pi ... % angular
                .* repmat(dateNormalised', 1, length(freqSpec));
            GSpec1 = [G1(:, 1:end - 2 * nOmega), ...
                          cos(phaseSpec), sin(phaseSpec)];
            GSpec2 = [G2(:, 1:end - 2 * nOmega), ...
                          cos(phaseSpec), sin(phaseSpec)];
            GSpec3 = [G3(:, 1:end - 2 * nOmega), ...
                          cos(phaseSpec), sin(phaseSpec)];
        else
            freqSpec = [];
            phaseSpec = [];
            GSpec1 = [];
            GSpec2 = [];
            GSpec3 = [];
        end

    end

    %% Compute the fits
    isSpecialterm = 1:nSleptCoeff == specialterms{1};
    % Since each Slepian coefficient has different errors, each will have a
    % different weighting matrix.  Thus we loop over the coefficients.
    parfor iSleptCoeff = 1:nSleptCoeff
        [sleptCoeffsSignal(:, iSleptCoeff), sleptCoeffsResid(:, iSleptCoeff), ...
             ftests(iSleptCoeff, :), extravalues(:, iSleptCoeff)] ...
            = slept2resid_fitslept(sleptCoeffs(:, iSleptCoeff), fitwhat, givenerrors(:, iSleptCoeff), ...
            isSpecialterm(iSleptCoeff), freq, phase, freqSpec, phaseSpec, dateNormalised, dateExtraNormalised, moredates, ...
            G1, G2, G3, GSpec1, GSpec2, GSpec3, nMonth);

    end

    %% Returning requested output
    varargout = {sleptCoeffsSignal, sleptCoeffsResid, ftests, extravalues};

    % Total combined fitting

    % If we have the parameters for this localization, and we requested the
    % total fit, then let's do that.

    if ~(nargout >= 5 && exist('CC', 'var') && exist('TH', 'var'))
        return
    end

    % Get the residual covariance
    [Cab] = slepresid2cov(sleptCoeffsResid);

    % Calculate the bandwdith for this basis
    L = sqrt(size(sleptCoeffs, 2)) - 1;
    % This should be an integer
    if (floor(L) ~= L)
        error('Something fishy about your L');
    end

    if iscell(TH)
        % Something like {'greenland' 0.5}
        if length(TH) >= 3
            XY = feval(TH{:});
        else
            XY = feval(TH{1}, 0, TH{2}, moreRegionSpecs{:});
        end

    else
        % Coordinates or a string, either works
        XY = TH;
    end

    % Calculate the Shannon number for this basis
    defval('N', round((L + 1) ^ 2 * spharea(XY)));

    % Make the coefficients with reference to some mean
    % If they already are, then this won't matter
    sleptdelta = sleptCoeffs(1:nMonth, :) ...
        - repmat(mean(sleptCoeffs(1:nMonth, :), 1), nMonth, 1);

    % COMBINE

    % We want to take the Slepian functions and combine them to get total mass.
    % For signal, this means integrating the functions and adding them.  For
    % the error, this means using the error propogation equation, where we
    % compute (int)*(covar)*(int)'.  Since the slepcoffs are constants that
    % just come forward, we can do the integration of the eigenfunctions
    % first (and once for each function), then multiply by slepcoffs to
    % get the monthly values.  This is much faster.

    if iscell(TH) && length(TH) > 2
        eigfunINT = integratebasis_new(CC, TH, N);
    else
        eigfunINT = integratebasis_new(CC, TH, N, "MoreRegionSpecs", moreRegionSpecs);
    end

    % Since Int should have units of (fn * m^2), need to go from fractional
    % sphere area to real area.  If the fn is surface density, this output is
    % in kilograms.  Then change the units from kg to Gt in METRIC tons
    eigfunINT = eigfunINT * 4 * pi * 6370000 ^ 2/10 ^ 3/10 ^ 9;
    functionintegrals = eigfunINT;

    % Now multiply by the appropriate slepcoffs to get the months
    % This becomes alpha by months
    %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
    %functimeseries = sleptdelta(:,1:N)';

    % Here do the total sum of the data
    eigfunINT = eigfunINT(:);
    eigfunINT = eigfunINT(1:N);
    total = eigfunINT' * sleptdelta(:, 1:N)';

    % Get the error
    thevars = diag(Cab(1:N, 1:N))';
    alphavar = eigfunINT .^ 2 .* thevars;
    % Now the combined error with covariance
    alphavarall = eigfunINT' * Cab(1:N, 1:N) * eigfunINT;

    % FITTING

    % We have uniform estimated error, which will be different than the polyfit
    % estimated residuals because ours account for sinusoidal signals.  So
    % pass the new error to our function for replacement, so
    % that the fitting confidence intervals reflect that

    [fit, delta, totalparams, paramerrors] = ...
        timeseriesfit([dates', total'], alphavarall, 1, 1);

    % Make a matrix for the line, and 95% confidence in the fit
    totalfit = [dates', fit, delta];

    % Make the error valid for a year
    totalparamerrors = paramerrors * 365;

    if strcmp(OutputUnit, 'year')
        totalparams = totalparams * 365;
    end

    % Collect the expanded output
    varargout = ...
        {sleptCoeffsSignal, sleptCoeffsResid, ftests, extravalues, ...
         total, alphavarall, totalparams, totalparamerrors, totalfit, ...
         functionintegrals, alphavar};

end
