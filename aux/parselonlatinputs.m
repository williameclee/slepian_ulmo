function varargout = parselonlatinputs(varargin)
    % Find longitude and latitude
    if isnumeric(varargin{1})

        if isvector(varargin{1}) && ...
                length(varargin) >= 2 && isvector(varargin{2})
            % First input is lon, second is lat
            lon = varargin{1}(:);
            lat = varargin{2}(:);

            % Check sizes
            if ~isequal(size(lon), size(lat))
                error('Longitude and latitude must have the same size.')
            end

            % Remove lonlat from varargin
            varargin(1:2) = [];

        elseif any(size(varargin{1}) == 2)
            % First input is lonlat
            lonlat = varargin{1};

            if size(lonlat, 1) == 2
                lonlat = lonlat';
            end

            lon = lonlat(:, 1);
            lat = lonlat(:, 2);

            varargin(1) = [];

        else
            error('Invalid input argument')
        end

    elseif isa(varargin{1}, 'polyshape')
        % First input is a polyshape
        [lon, lat] = poly2xy(varargin{1});
        varargin(1) = [];
    elseif isa(varargin{1}, 'GeoDomain')
        % First input is a GeoDomain
        domain = varargin{1};
        lonlat = domain.Lonlat;
        lon = lonlat(:, 1);
        lat = lonlat(:, 2);
        varargin(1) = [];
    else
        error('Invalid input argument type for the first argument')
    end

    varargout = {lon, lat, varargin};
end
