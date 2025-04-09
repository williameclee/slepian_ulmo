%% REMOVEPOLARCAPS
% Truncate a polygon to between the specified latitude band.
%
% Input arguments
%   poly - The polygon (POLYSHAPE) to truncate
%   inclang - The inclination angle of the polar caps
%       The inclination angle can be specified as a scalar or a 1-by-2
%       vector.
%       - If a scalar, the inclination angle is applied to both the north
%           and south polar caps.
%       - If a 1-by-2 vector, the first element is the inclination angle of
%           the north polar cap, and the second element is the inclination
%           angle of the south polar cap.
%
% Last modified by
%   2025/04/09, williameclee@arizona.edu (@williameclee)

function poly = removepolarcaps(poly, inclang, ~)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');

    if isscalar(inclang)
        inclang = [-inclang, inclang];
    end

    %% Generating the polar caps
    eps = 3;
    minLon = min(poly.Vertices(:, 1)) - eps;
    maxLon = max(poly.Vertices(:, 1)) + eps;

    if inclang(1) > -90
        [XcapS, YcapS] = poly2ccw( ...
            [minLon; maxLon; maxLon; minLon], ...
            [-90 - eps; -90 - eps; inclang(1); inclang(1)]);
        capS = polyshape(XcapS, YcapS);
    else
        capS = polyshape();
    end

    if inclang(2) < 90
        [XcapN, YcapN] = poly2ccw( ...
            [minLon; maxLon; maxLon; minLon], ...
            [inclang(2); inclang(2); 90 + eps; 90 + eps]);
        capN = polyshape(XcapN, YcapN);
    else
        capN = polyshape();
    end

    %% Removing the polar caps and returning the result
    poly = subtract(poly, union(capN, capS));
end
