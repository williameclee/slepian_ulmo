function [upscale, inclang, buf, moreBuf, forceNew, lonOrigin] = ...
        parsecoastinputs(inputArguments, varargin)
    %% Assigning default values
    d = inputParser;
    addRequired(d, 'inputArguments');
    addOptional(d, 'DefaultUpscale', 0);
    addOptional(d, 'DefaultBuffer', 0);
    addOptional(d, 'DefaultInclang', 90);
    addOptional(d, 'DefaultMoreBuffer', []);
    addOptional(d, 'DefaultForceNew', false);
    addOptional(d, 'DefaultLongitudeOrigin', 0);
    parse(d, inputArguments, varargin{:});
    defaultUpscale = d.Results.DefaultUpscale;
    defaultInclang = d.Results.DefaultInclang;
    defaultBuffer = d.Results.DefaultBuffer;
    defaultMoreBuffer = d.Results.DefaultMoreBuffer;
    defaultForceNew = d.Results.DefaultForceNew;
    defaultLongitudeOrigin = d.Results.DefaultLongitudeOrigin;

    %% Parsing the 'real' inputs
    p = inputParser;
    p.KeepUnmatched = true;
    addOptional(p, 'Upscale', defaultUpscale, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', defaultBuffer, @isnumeric);
    addOptional(p, 'Inclang', defaultInclang, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'MoreBuffer', defaultMoreBuffer, @iscell);
    addParameter(p, 'ForceNew', defaultForceNew, ...
        @(x) islogical(x) || x == 1 || x == 0);
    addParameter(p, 'LongitudeOrigin', defaultLongitudeOrigin);
    parse(p, inputArguments{:});

    %% Assigning the parsed values
    upscale = p.Results.Upscale;

    if isempty(upscale)
        upscale = defaultUpscale;
    end

    inclang = p.Results.Inclang;

    if any(isempty(inclang)) || any(isnan(inclang))
        inclang = defaultInclang;
    end

    buf = p.Results.Buffer;

    if isempty(buf)
        buf = defaultBuffer;
    end

    moreBuf = p.Results.MoreBuffer;

    if isempty(moreBuf)
        moreBuf = defaultMoreBuffer;
    end

    forceNew = logical(p.Results.ForceNew);
    lonOrigin = p.Results.LongitudeOrigin;
end
