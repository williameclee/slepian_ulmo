%% REMOVEDUPLICATEVERTICES
% Removes duplicate vertices from a polyline.
%
% Syntax
%   XYs = removeduplicatevertices(XY)
%   XYs = removeduplicatevertices(X, Y)
%   [Xs, Ys] = removeduplicatevertices(__)
%
% Input arguments
%   XY - The vertices of the polyline
%       The vertices are stored in a 2-column matrix or a polyshape.
%   X, Y - The x- and y-coordinates of the vertices
%
% Output arguments
%   XYs - The vertices of the polyline without duplicates
%   Xs, Ys - The x- and y-coordinates of the vertices without duplicates
%
% Last modified by
%   2024/08/12, williameclee@arizona.edu (@williameclee)

function varargout = removeduplicatevertices(varargin)

    if nargin == 0 || nargin > 2
        error('Invalid number of inputs.')
    elseif nargin == 1

        if isa(varargin{1}, 'polyshape')
            varargin{1} = varargin{1}.Vertices;
        end

        XY = varargin{1};

        if ~any(size(XY) == 2)
            error('The input must be a 2-column matrix.')
        elseif size(XY, 1) == 2
            XY = XY';
        end

    elseif nargin == 2
        X = varargin{1}(:);
        Y = varargin{2}(:);

        if ~isequal(size(X), size(Y))
            error('X and Y must have the same size.')
        end

        XY = [X, Y];
    end

    isDuplicate = ...
        [false; all(XY(1:end - 1, :) == XY(2:end, :), 2)];
    XYs = XY(~isDuplicate, :);

    if nargout == 0
        figName = sprtinf('Vertices (%s)', upper(mfilename));
        figure(999)
        set(gcf, "Name", figName, "NumberTitle", 'off')

        plot(XYs(:, 1), XYs(:, 2), 'k.-')

        varargout = {XYs};
        return
    end

    if nargin == 1
        varargout = {XYs};
    elseif nargin == 2
        varargout = {XYs(:, 1), XYs(:, 2)};
    else
        error('Invalid number of outputs.')
    end

end
