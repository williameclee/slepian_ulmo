classdef GeoDomain

    properties
        Domain string;
        Upscale(1, 1) double {mustBeNumeric};
        Buffer(1, 1) double {mustBeNumeric};
        Latlim(1, 2) double {mustBeNumeric};
        MoreBuffers cell;
        NearBy(1, 1) logical;
    end

    methods

        % Constructor
        function obj = GeoDomain(varargin)
            % Parse inputs
            Inputs = parseinputs(varargin{:});

            % Empty domain (for preallocation)
            if isempty(Inputs.Domain)
                return
            end

            % Find out how many elements are should be in the object
            % Also turn the inputs into the right size
            [nObjects, domain, upscale, buf, latlim, moreBuf, nearby] = ...
                noutputs(lower(Inputs.Domain), Inputs.Upscale, ...
                Inputs.Buffer, Inputs.Latlim, ...
                Inputs.MoreBuffers, Inputs.NearBy);

            % Whether to use default parameters
            useDefaultParams = Inputs.DefaultParams;

            % Preallocation
            for i = 1:nObjects
                obj(i, 1) = GeoDomain(''); %#ok<AGROW>
            end

            obj(nObjects) = GeoDomain('');

            % Assign the domain values for each object
            for i = 1:nObjects
                obj(i) = assigndomainvalues(domain(i), upscale(i), ...
                    buf(i), latlim(i, :), nearby(i), moreBuf{i}, ...
                    useDefaultParams);
            end

        end

        % The identifier for file names, etc.
        function id = Id(obj)
            id = ...
                [capitalise(char(obj.Domain)), dataattrchar( ...
                 "Upscale", obj.Upscale, "Buffer", obj.Buffer, ...
                 "Latlim", obj.Latlim, "MoreBuffers", obj.MoreBuffers)];
        end

        % Display name
        function displayName = DisplayName(obj, varargin)
            p = inputParser;
            addOptional(p, 'Format', 'short', ...
                @(x) (ischar(x) || isstring(x)) && ...
                ismember(lower(x), {'short', 'long'}));
            parse(p, varargin{:});
            format = p.Results.Format;

            displayName = domainname(obj.Domain, format);

        end

        % Longitude and latitude arrays
        function lonlat = Lonlat(obj, varargin)
            p = inputParser;
            addParameter(p, 'RotateBack', false, ...
                @(x) islogical(x) || isnumeric(x));
            parse(p, varargin{:});
            rotateBack = logical(p.Results.RotateBack);

            if ~strcmp(obj.Domain, 'antarctica') && rotateBack
                rotateBack = false;
                warning('GeoDomain:unsupportedRotateBack', ...
                'RotateBack is only supported for Antarctica');
            end

            lonlat = feval(obj.Domain, "Upscale", obj.Upscale, ...
                "Buffer", obj.Buffer, "Latlim", obj.Latlim, ...
                "MoreBuffers", obj.MoreBuffers, "RotateBack", rotateBack);

            if nargout > 0
                return
            end

            plotlonlat(lonlat, obj, rotateBack)

        end

        % Longitude only
        function lon = Lon(obj, varargin)
            lonlat = obj.Lonlat(varargin{:});
            lon = lonlat(:, 1);
        end

        % Latitude only
        function lat = Lat(obj, varargin)
            lonlat = obj.Lonlat(varargin{:});
            lat = lonlat(:, 2);
        end

        % Fractional area relative to the sphere
        function area = SphArea(obj)
            area = spharea(obj.Lonlat);
        end

        % Total area in square metres
        function area = Area(obj)
            earthRadius = 6371e3;
            area = spharea(obj.Lonlat) * (4 * pi * earthRadius ^ 2);
        end

    end

end

%% Subfunctions
% Parse inputs
function Inputs = parseinputs(varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p, 'Domain', ...
        @(x) ischar(x) || isstring(x) || iscell(x));
    addOptional(p, 'Upscale', [], @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', [], @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Latlim', [], ...
        @(x) (isnumeric(x) && length(x) <= 2) ...
        || isempty(x) || strcmpi(x, 'default'));
    addOptional(p, 'MoreBuffers', {}, ...
        @(x) iscell(x) || isempty(x) || strcmpi(x, 'default'));
    addOptional(p, 'NearBy', false, ...
        @(x) (isnumeric(x) && x <= 1) || islogical(x));
    addParameter(p, 'DefaultParams', false, @islogical);
    parse(p, varargin{:});

    Inputs = p.Results;
end

% Find out how many elements are should be in the object
function [nObjects, domain, upscale, buf, latlim, moreBuf, nearby] = ...
        noutputs(domain, upscale, buf, latlim, moreBuf, nearby)

    if iscell(domain) || isstring(domain)
        nDomains = length(domain);

        if iscell(domain)
            domain = string(domain(:));
        end

    else
        nDomains = 1;
        domain = string(domain);
    end

    nUpscales = length(upscale);

    if isempty(upscale)
        upscale = nan;
    elseif all(strcmpi(upscale, 'default'))
        upscale = nan;
    end

    nBuffers = length(buf);

    if isempty(buf)
        buf = nan;
    elseif all(strcmpi(buf, 'default'))
        buf = nan;
    end

    if isempty(latlim) || strcmpi(latlim, 'default')
        latlim = [nan, nan];
        nLatlims = 1;
    else

        if ~isscalar(latlim)

            if size(latlim, 1) > 1
                nLatlims = size(latlim, 1);

                if size(latlim, 2) == 1
                    latlim = latlim * [-1, 1];
                end

            elseif size(latlim, 1) == 2 && size(latlim, 2) > 2
                nLatlims = size(latlim, 2);
                latlim = latlim';
            elseif isequal(size(latlim), [1, 2])
                nLatlims = 1;
            else
                error('The latitude limits must be a 1x2 or 2xN array');
            end

        else
            nLatlims = 1;
            latlim = [-1, 1] * latlim;
        end

    end

    if isempty(moreBuf)
        nMoreBuffers = 1;
        moreBuf = {''};
    elseif strcmpi(moreBuf, 'default')
        moreBuf = "default";
    else

        if iscell(moreBuf{1})
            nMoreBuffers = length(moreBuf);
        else
            nMoreBuffers = 1;
        end

    end

    nNearby = length(nearby);

    nObjects = unique([nDomains, nUpscales, nBuffers, nLatlims, ...
                           nMoreBuffers, nNearby]);
    nObjects = nObjects(nObjects >= 1);

    if length(nObjects) > 2
        error('The number of elements in the input arguments must be the same');
    else
        nObjects = nObjects(end);
    end

    if nObjects == 1
        moreBuf = {moreBuf};
        return
    end

    if nDomains < nObjects
        domain = repmat(domain, nObjects, 1);
    end

    if nUpscales < nObjects
        upscale = repmat(upscale, nObjects, 1);
    end

    if nBuffers < nObjects
        buf = repmat(buf, nObjects, 1);
    end

    if nLatlims < nObjects
        latlim = repmat(latlim, nObjects, 1);
    end

    if nMoreBuffers < nObjects
        moreBuf = repmat(moreBuf, nObjects, 1);
    end

    if nNearby < nObjects
        nearby = repmat(nearby, nObjects, 1);
    end

end

% Assign the domain values (in case of multiple domains)
function obj = ...
        assigndomainvalues(domain, upscale, buf, latlim, nearby, ...
        moreBuf, useDefaultParams)
    defaultUpscale = 0;
    dafaultBuffer = 0;
    defaultLatlim = [-90, 90];
    defaultNearby = false;
    obj = GeoDomain('');
    obj.Domain = domain;

    if exist(obj.Domain, 'file') ~= 2
        warning('GeoDomain:domainNotFound', ...
            'Domain %s might not exist', upper(obj.Domain));
    end

    if isnan(upscale)
        upscale = [];
    end

    obj.Upscale = conddefval(upscale, defaultUpscale);

    if obj.Upscale < 0
        warning('GeoDomain:invalidUpscaleValue', ...
            ['Upscale must be non-negative', newline, ...
             sprintf('Changed from %d to 0', obj.Upscale)]);
        obj.Upscale = 0;
    elseif obj.Upscale == 1
        obj.Upscale = 0;
    end

    if isnan(buf)
        buf = [];
    end

    obj.Buffer = conddefval(buf, dafaultBuffer);

    if strcmpi(latlim, 'default') || any(isnan(latlim)) || ...
            (useDefaultParams && isempty(p.Results.Latlim))

        switch obj.Domain
            case {'spacific', 'satlantic', 'indian'}
                latlim = [-60, 60];
            otherwise
                latlim = [-90, 90];
        end

    else
        latlim = conddefval(latlim, defaultLatlim);
    end

    if latlim(1) == latlim(2)
        error('%s\nThe given range is: [%i, %i]', ...
            'The upper and lower bounds of the latitude limits must be different', ...
            latlim(1), latlim(2));
    else
        obj.Latlim = latlim;
    end

    if any(strcmpi(moreBuf, 'default')) || ...
            (useDefaultParams && (isempty(moreBuf) || ...
            all(cellfun(@isempty, moreBuf))))

        switch obj.Domain
            case {'oceans', 'npacific', 'indian'}
                moreBuf = {'earthquakes', 10};
            otherwise
                moreBuf = {''};
        end

    end

    % Sort the more buffers
    if ~all(cellfun(@isempty, moreBuf))
        moreBufDomain = lower(moreBuf(1:2:end));
        moreBufWidth = moreBuf(2:2:end);
        [moreBuf(1:2:end), moreBufSortId] = sort(moreBufDomain);
        moreBuf(2:2:end) = moreBufWidth(moreBufSortId);
    else
        moreBuf = {};
    end

    obj.MoreBuffers = moreBuf;

    obj.NearBy = conddefval(nearby, defaultNearby);

    if ~strcmp(obj.Domain, 'greenland') && obj.NearBy
        obj.NearBy = false;
        warning('GeoDomain:unsupportedNearBy', ...
            ['Option NEARBY is only supported for Greenland', ...
             newline, 'Changed from true to false']);
    end

end

% Plot a dirty map of the boundary
function plotlonlat(lonlat, obj, rotateBack)
    figure(999)
    set(gcf, 'Name', obj.DisplayName, 'NumberTitle', 'off')
    clf

    dirtymap(lonlat)

    if obj.Buffer == 0
        return
    end

    lonlatnb = feval(obj.Domain, "Upscale", obj.Upscale, "Buffer", 0, ...
        "Inclang", obj.Latlim, "MoreBuffer", obj.MoreBuffers, ...
        "RotateBack", rotateBack);
    hold on
    plot(lonlatnb(:, 1), lonlatnb(:, 2), 'b', ...
        "DisplayName", "Unbuffered")
    plot(lonlat(:, 1), lonlat(:, 2), 'k', ...
        "DisplayName", "Buffered")
    hold off
end
