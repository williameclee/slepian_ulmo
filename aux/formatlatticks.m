function varargout = formatlatticks(varargin)

    if nargin == 0
        lat = yticks;
    elseif nargin == 1
        lat = varargin{1};
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

    if nargout == 0
        yticklabels(latL)
        return
    else
        varargout = {latL};
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
