function [coastXY, coastPoly] = buffer4oceans(coast, varargin)
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    % Parse inputs
    p = inputParser;
    addRequired(p, 'Coast');
    addOptional(p, 'MoreBuffers', []);
    addOptional(p, 'LonOrigin', 0);
    parse(p, coast, varargin{:});
    coast = p.Results.Coast;
    moreBufs = p.Results.MoreBuffers;
    lonOrigin = p.Results.LonOrigin;

    %% Buffer continents with different buffers
    if isa(coast, 'polyshape')
        coastPoly = coast;
    else
        coastPoly = polyshape(coast);
    end

    if isempty(moreBufs)
        coastXY = closecoastline(coastPoly.Vertices);
        return
    end

    for i = 1:length(moreBufs) / 2
        moreCoastXY = feval(moreBufs{i * 2 - 1}, [], moreBufs{i * 2});
        [moreCoastY, moreCoastX] = flatearthpoly( ...
            moreCoastXY(:, 2), moreCoastXY(:, 1), lonOrigin);
        moreCoastPoly = polyshape(moreCoastX, moreCoastY);
        coastPoly = union(coastPoly, moreCoastPoly);
    end

    coastXY = closecoastline(coastPoly.Vertices);

end
