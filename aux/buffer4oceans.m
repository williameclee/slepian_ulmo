function coastPoly = buffer4oceans(coast, buf, varargin)
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    % Parse inputs
    p = inputParser;
    addRequired(p, 'Coast');
    addOptional(p, 'buf', []);
    addOptional(p, 'MoreBuffer', []);
    addOptional(p, 'LongitudeOrigin', 0);
    parse(p, coast, buf, varargin{:});
    coast = p.Results.Coast;
    % buf = p.Results.buf;
    moreBuf = p.Results.MoreBuffer;
    lonOrigin = p.Results.LongitudeOrigin;

    %% Buffer continents with different buffers
    if isa(coast, 'polyshape')
        coastPoly = coast;
    else
        coastPoly = polyshape(coast);
    end

    if isempty(moreBuf)
        return
    end

    for i = 1:length(moreBuf) / 2
        moreCoastXY = feval(moreBuf{i * 2 - 1}, [], moreBuf{i * 2});
        [moreCoastY, moreCoastX] = flatearthpoly( ...
            moreCoastXY(:, 2), moreCoastXY(:, 1), lonOrigin);
        moreCoastPoly = polyshape(moreCoastX, moreCoastY);
        coastPoly = union(coastPoly, moreCoastPoly);
    end

end
