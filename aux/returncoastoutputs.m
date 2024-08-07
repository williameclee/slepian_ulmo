function vout = returncoastoutputs(nArgout, XY, varargin)
    p = inputParser;
    addRequired(p, 'nArgout', @isnumeric);
    addRequired(p, 'XY', @isnumeric);
    addOptional(p, 'MapTitle', '', @ischar);
    addParameter(p, 'FigureNumber', 999, @isnumeric);
    parse(p, nArgout, XY, varargin{:});

    nArgout = p.Results.nArgout;
    XY = p.Results.XY;
    mapTitle = p.Results.MapTitle;
    figNum = p.Results.FigureNumber;

    vout = {XY};

    if nArgout > 0
        return
    end

    figure(figNum)
    clf
    set(gcf, 'Name', mapTitle, 'NumberTitle', 'off')

    dirtymap(XY)
    hold on
    plot(XY(:, 1), XY(:, 2), 'k.-')
    hold off
    title(['Point count: ', num2str(size(XY, 1))])

end
