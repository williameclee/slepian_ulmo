function varargout = formatlonticks(varargin)

    if nargin == 0
        lon = xticks;
    elseif nargin == 1
        lon = varargin{1};
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

    if nargout == 0
        xticklabels(lonL)
        return
    else
        varargout = {lonL};
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
