function Plm = sphgaussfilt(Plm, varargin)
    ip = inputParser;
    addRequired(ip, 'Plm', @(x) (isnumeric(x) && size(x, 2) == 4));
    addOptional(ip, 'Sigma', 5, @isnumeric);
    addOptional(ip, 'Unit', 'degree', @(x) ischar(validatestring(x, {'degree', 'radian'})));
    parse(ip, Plm, varargin{:});
    Plm = ip.Results.Plm;
    sigma = ip.Results.Sigma;
    unit = ip.Results.Unit;

    switch unit
        case 'degree'
            sigma = 180 / sqrt(2) / sigma;
    end

    if ismatrix(Plm)
        l = Plm(:, 1);
        Plm(:, 3:4) = Plm(:, 3:4) .* exp(-l .* (l + 1) / sigma ^ 2/2);
    elseif ndims(Plm) == 3
        l = Plm(:, 1, :);
        Plm(:, 3:4, :) = Plm(:, 3:4, :) .* exp(-l .* (l + 1) / sigma ^ 2/2);
    else
        error('Invalid input dimensions');
    end

end