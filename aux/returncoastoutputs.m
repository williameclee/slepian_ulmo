function vout = returncoastoutputs(nOut, XY, oceanPoly, varargin)
    p = inputParser;
    addRequired(p, 'nOut', @isnumeric);
    addRequired(p, 'XY', @isnumeric);
    addRequired(p, 'Polygon', @(x) isa(x, 'polyshape'));
    addOptional(p, 'CallerName', '', @ischar);
    addParameter(p, 'FigureNumber', 999, @isnumeric);
    parse(p, nOut, XY, oceanPoly, varargin{:});

    nOut = p.Results.nOut;
    XY = p.Results.XY;
    oceanPoly = p.Results.Polygon;
    callerName = p.Results.CallerName;
    figNum = p.Results.FigureNumber;

    vout = {XY, oceanPoly};

    if nOut > 0
        return
    end

    %% Plotting the boundary
    if isempty(callerName)
        st = dbstack;
        callerName = st(2).name;
    end

    figTitle = sprintf('%s (%s)', domainname(callerName, 'long'), upper(callerName));

    figure(figNum)
    clf
    set(gcf, 'Name', figTitle, 'NumberTitle', 'off')

    plotqdm(XY, 'k')
    % hold on
    % plot(XY(:, 1), XY(:, 2), 'k.-')
    % hold off
    title(['Point count: ', num2str(size(XY, 1))])

end
