function poly = removepolarcaps(poly, inclang, lonOrigin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse inputs
    p = inputParser;
    addRequired(p, 'poly');
    addOptional(p, 'inclang', 90);
    addOptional(p, 'lonOrigin', 0);

    if isscalar(inclang)
        inclang = [-inclang, inclang];
    end

    %% Generating the polar caps
    if inclang(1) > -90
        [XcapS, YcapS] = poly2ccw( ...
            [lonOrigin - 180 - 10; lonOrigin + 180 + 10; lonOrigin + 180 + 10; lonOrigin - 180 - 10], ...
            [-90 - 10; -90 - 10; inclang(1); inclang(1)]);
        capS = polyshape(XcapS, YcapS);
    else
        capS = polyshape();
    end

    if inclang(2) < 90
        [XcapN, YcapN] = poly2ccw( ...
            [lonOrigin - 180 - 10; lonOrigin + 180 + 10; lonOrigin + 180 + 10; lonOrigin - 180 - 10], ...
            [inclang(2); inclang(2); 90 + 10; 90 + 10]);
        capN = polyshape(XcapN, YcapN);
    else
        capN = polyshape();
    end

    %% Removing the polar caps and returning the result
    poly = subtract(poly, union(capN, capS));

end
