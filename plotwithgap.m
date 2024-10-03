function plotwithgap(x, y, tol, varargin)
    [x, y] = splitbygap(x, tol, y);
    plot(x, y, varargin{:});
end
