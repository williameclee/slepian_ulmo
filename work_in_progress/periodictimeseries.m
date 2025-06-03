%% PERIODICTIMESERIES
%
% Authored by:
%   2025/06/03, williameclee@arizona.edu (@williameclee)

function [polyCoeffs, periodicCoeffs, dataFit] = periodictimeseries(t, varargin)
    ip = inputParser;
    ip.addRequired('t', ...
        @(x) (isnumeric(x) && (isvector(x) || ismatrix(x))) || ...
        ((isdatetime(x) || isduration(x)) && isvector(x)));
    ip.addOptional('x', [], @(x) isnumeric(x) && isvector(x));
    ip.addOptional('p', 2, @(x) isnumeric(x) && isscalar(x));
    ip.addOptional('periods', [], @(x) isnumeric(x) || isduration(x));
    ip.addParameter('PeriodicFormat', 'cos-sin', @(x) ischar(validatestring(x, {'sin-cos', 'cos-sin', 'amp-phase'})));
    ip.parse(t, varargin{:});

    t = ip.Results.t;
    x = ip.Results.x;
    p = ip.Results.p;
    periods = ip.Results.periods;
    periodicFormat = ip.Results.PeriodicFormat;

    if isempty(x) && isnumeric(t) && any(size(t) == 2)

        if size(t, 2) == 2
            x = t(:, 2);
            t = t(:, 1);
        elseif size(t, 1) == 2
            x = t(2, :);
            t = t(1, :);
        else
            error('Invalid input dimensions');
        end

    elseif isempty(x)
        error('No data provided');
    end

    N = numel(t);

    if N ~= numel(x)
        error('dates and data must have the same number of elements');
    end

    t = t(:);
    x = x(:);

    isTime = false;

    if isdatetime(t)
        isTime = true;
        t = years(t - datetime(year(t(1)), 1, 1));
    end

    if ~isempty(periods) && isduration(periods)
        periods = years(periods);
    end

    X = zeros([N, p + 1 + 2 * numel(periods)]);
    X(:, 1:p + 1) = t .^ (0:p);

    for i = 1:numel(periods)

        switch periodicFormat
            case "sin-cos"
                X(:, p + 1 + (i - 1) * 2 + [1, 2]) = ...
                    [sin(2 * pi * t / periods(i)), cos(2 * pi * t / periods(i))];
            otherwise
                X(:, p + 1 + (i - 1) * 2 + [1, 2]) = ...
                    [cos(2 * pi * t / periods(i)), sin(2 * pi * t / periods(i))];
        end

    end

    coeffs = X \ x;
    dataFit = X * coeffs;

    polyCoeffs = coeffs(1:p + 1);
    periodicCoeffs = coeffs(p + 2:end);

    if strcmp(periodicFormat, "amp-phase")
        periodicCoeffs = reshape(periodicCoeffs, [2, numel(periods)]);
        periodicCoeffs = ...
            [sqrt(sum(periodicCoeffs .^ 2, 1)); ...
             wrapTo2Pi(atan2(periodicCoeffs(2, :), periodicCoeffs(1, :)))];
        periodicCoeffs = periodicCoeffs';

        if isTime
            periodicCoeffs(:, 2) = periodicCoeffs(:, 2) / (2 * pi) * days(years(1));
        end

    end

end
