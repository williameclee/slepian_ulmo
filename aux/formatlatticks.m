%% FORMATLATTICKS
% Format latitude tick labels for a plot.
%
% Syntax
%   formatlatticks
%       Formats the y-axis tick labels of the current plot as latitude
%       ticks.
%   formatlatticks(lat)
%       Formats the y-axis tick labels of the current plot as latitude
%       ticks with the specified values.
%   formatlatticks(ax)
%       Formats the y-axis tick labels of the specified axes as latitude
%       ticks.
%   latL = formatlatticks(__)
%       Returns the formatted latitude tick labels as a cell array of
%       character vectors.
%
% Input arguments
%   lat - Latitude ticks in degrees
%   ax - Axes object to format
%
% Output arguments
%   latL - Formatted latitude tick labels
%
% Last modified by
%   2024/08/14, williameclee@arizona.edu (@williameclee)

function varargout = formatlatticks(varargin)

    if nargin == 0
        lat = yticks;
    elseif nargin == 1

        if isnumeric(varargin{1})
            lat = varargin{1};
        elseif isa(varargin{1}, 'matlab.graphics.axis.Axes')
            lat = yticks(varargin{1});
        else
            error('Invalid input type %s', class(varargin{1}));
        end

    else
        error('Invalid number of arguments')
    end

    if isscalar(lat)
        latL = formatlontick(lat);
    elseif isnumeric(lat)
        latL = arrayfun(@formatlattick, lat, 'UniformOutput', false);
    elseif iscell(lat)
        latL = cellfun(@formatlattick, lat, 'UniformOutput', false);
    else
        error('Invalid input type %s', class(lat));
    end

    if nargout > 0
        varargout = {latL};
        return
    end

    clear varargout

    if ~isempty(varargin) && ...
            isa(varargin{1}, 'matlab.graphics.axis.Axes')
        yticklabels(varargin{1}, latL)
    else
        yticklabels(latL)
    end

end

%% Subfunctions
function latL = formatlattick(lat)

    if lat > 0
        latL = sprintf('%d°N', lat);
    elseif lat < 0
        latL = sprintf('%d°S', -lat);
    else
        latL = '0°';
    end

end
