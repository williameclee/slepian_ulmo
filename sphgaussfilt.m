%% SPHGAUSSFILT
% Applys a Gaussian filter to a spherical harmonic field.
%
% Syntax
%   Plm = sphgaussfilt(Plm, sigma)
%   Plm = sphgaussfilt(Plm, sigma, unit)
%
% Input arguments
%   Plm - Spherical harmonic coefficients
%   sigma - Standard deviation of the Gaussian filter
%   unit - Unit of the standard deviation
%       'degree' - The standard deviation is in degrees
%       'radian' - The standard deviation is in radians
%       'km' - The standard deviation is in kilometres
%       The default unit is 'degree'
%
% Output arguments
%   Plm - Filtered spherical harmonic coefficients
%
% Last modified by
%   2025/03/18, williameclee@arizona.edu (@williameclee)
%   2024/11/20, williameclee@arizona.edu (@williameclee)

function Plm = sphgaussfilt(Plm, varargin)
    ip = inputParser;
    addRequired(ip, 'Plm', @(x) (isnumeric(x) && size(x, 2) == 4));
    addOptional(ip, 'Sigma', 5, @isnumeric);
    addOptional(ip, 'Unit', 'degree', @(x) ischar(validatestring(x, {'degree', 'radian', 'km'})));
    parse(ip, Plm, varargin{:});
    Plm = ip.Results.Plm;
    sigma = ip.Results.Sigma;
    unit = ip.Results.Unit;

    if ismatrix(Plm)
        l = Plm(:, 1);
    elseif ndims(Plm) == 3
        l = Plm(:, 1, :);
    end

    switch unit
        case 'degree'
            sigma = 180 / sqrt(2) / sigma;
            gfilter = exp(-l .* (l + 1) / sigma ^ 2/2);
        case 'radian'
            gfilter = exp(-l .* (l + 1) / sigma ^ 2/2);
        case 'km' % modified from Adhikari et al. (2019)
            % sigma = sigma/6371;
            gfilter = exp(- l .* (l + 1) * (sigma / 6371) ^ 2 / (4 * log(2)));
    end

    if ismatrix(Plm)
        Plm(:, end - 1:end) = Plm(:, end - 1:end) .* gfilter;
    elseif ndims(Plm) == 3
        Plm(:, end - 1:end, :) = Plm(:, end - 1:end, :) .* gfilter;
    else
        error('Invalid input dimensions');
    end

end
