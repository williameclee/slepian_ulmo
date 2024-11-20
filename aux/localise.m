%% Localise
% Localise spherical harmonic coefficients to a domain.
% This goal is very similar to that of Slepian functions', but this function don't truncate anything.
%
% Syntax
%   Plm = localise(Plm, domain, L)
%   Plm = localise(Plm, domain, L, "Inverse", true)
%   Plm = localise(Plm, "K", K)
%   [Plm, K] = localise(Plm, __)
%
% Input arguments
%   Plm - Spherical harmonic coefficients
%       The coefficients are in the form lmcosi, but the first two columns (degree and order) can be omitted.
%       Can be three-dimensional, where the third dimension is the data.
%   domain - Geographic domain
%       A geographic domain (GeoDomain object).
%   L - Bandwidth of the window
%       The default value is the degree of the data.
%   "Inverse" - Invert the kernel
%       The default value is false.
%   "K" - Kernel matrix
%       The default value is computed using KERNELCP_NEW.
%       For large L, loading the kernel matrix can be time-consuming.
%
% Output arguments
%   Plm - Localised spherical harmonic coefficients
%   K - Kernel matrix
%
% See also
%   KERNELCP_NEW
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)

function [Plm, K] = localise(Plm, domain, L, varargin)
    ip = inputParser;
    addRequired(ip, 'Plm', ...
        @(x) isnumeric(x) && (size(x, 2) == 4 || size(x, 2) == 2));
    addRequired(ip, 'domain', @(x) isa(x, 'GeoDomain'));
    addRequired(ip, 'L', @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'Inverse', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'K', [], @ismatrix);
    parse(ip, Plm, domain, L, varargin{:});
    Plm = ip.Results.Plm;
    domain = ip.Results.domain;
    L = ip.Results.L;
    isInverted = ip.Results.Inverse;
    K = ip.Results.K;

    is3d = ndims(Plm) == 3;

    if is3d
        nData = size(Plm, 3);
    end

    includesDO = size(Plm, 2) == 4;

    if includesDO
    end

    j2 = kernelorder(L);

    if includesDO

        if is3d
            Plm = Plm(:, 3:4, :);
        else
            Plm = Plm(:, 3:4);
        end

    end

    if size(Plm, 1) < addmup(L)

        if is3d
            Plm(addmup(L), 2, 1) = 0;
        else
            Plm(addmup(L), 2) = 0;
        end

    elseif size(Plm, 1) > addmup(L)

        if is3d
            Plm = Plm(1:addmup(L), :, :);
        else
            Plm = Plm(1:addmup(L), :);
        end

    end

    if is3d
        Plm = reshape(Plm, [size(Plm, 1) * size(Plm, 2), size(Plm, 3)]);
        Plms = Plm(j2, :);
    else
        Plms = Plm(j2);
    end

    if isempty(K)
        K = kernelcp_new(L, domain, "BeQuiet", true);

        if isInverted
            K = eye(size(K)) - K;
        end

    end

    ofun = K * Plms;

    if is3d
        Plm = zeros([addmup(L) * 2, nData]);
        Plm(j2, :) = ofun;
        Plm = reshape(Plm, [addmup(L), 2, nData]);
    else
        Plm = zeros([addmup(L), 2]);
        Plm(j2) = ofun;
    end

    if includesDO

        [order, degree] = addmon(L);

        if is3d
            Plm = cat(2, repmat(degree, [1, 1, nData]), repmat(order, [1, 1, nData]), Plm);
        else
            Plm = [degree, order, Plm];
        end

    end

end
