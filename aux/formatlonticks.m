%% FORMATLONTICKS
% Format longitude tick labels for a plot.
%
% Syntax
%   formatlonticks
%       Formats the x-axis tick labels of the current plot as longitude
%       ticks.
%   formatlonticks(lon)
%       Formats the x-axis tick labels of the current plot as longitude
%       ticks with the specified values.
%   formatlonticks(ax)
%       Formats the x-axis tick labels of the specified axes as longitude
%       ticks.
%   lonL = formatlonticks(__)
%       Returns the formatted longitude tick labels as a cell array of
%       character vectors.
%
% Input arguments
%   lon - Longitude ticks in degrees
%   ax - Axes object to format
%
% Output arguments
%   lonL - Formatted longitude tick labels
%
% Last modified by
%   2024/08/14, williameclee@arizona.edu (@williameclee)

function varargout = formatlonticks(varargin)

    if nargin == 0
        lon = xticks;
    elseif nargin == 1

        if isnumeric(varargin{1})
            lon = varargin{1};
        elseif isa(varargin{1}, 'matlab.graphics.axis.Axes')
            lon = xticks(varargin{1});
        else
            error('Invalid input type %s', class(varargin{1}));
        end

    else
        error('Invalid number of arguments')
    end

    if isscalar(lon)
        lonL = formatlontick(lon);
    elseif isnumeric(lon)
        lonL = arrayfun(@formatlontick, lon, 'UniformOutput', false);
    elseif iscell(lon)
        lonL = cellfun(@formatlontick, lon, 'UniformOutput', false);
    else
        error('Invalid input type %s', class(lon));
    end

    if nargout > 0
        varargout = {lonL};
        return
    end

    clear varargout
    
    if ~isempty(varargin) && ...
            isa(varargin{1}, 'matlab.graphics.axis.Axes')
        xticklabels(varargin{1}, lonL)
    else
        xticklabels(lonL)
    end

end

%% Subfunctions
function lonL = formatlontick(lon)
    lon = mod(lon, 360);

    if lon < 180
        lonL = sprintf('%d°E', lon);
    elseif lon > 180
        lonL = sprintf('%d°W', 360 - lon);
    else
        lonL = '180°';
    end

end
