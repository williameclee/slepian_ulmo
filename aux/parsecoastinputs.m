%% PARSECOASTINPUTS
% Parses the inputs for the ocean domain functions.
%
% Syntax
%   parsecoastinputs(inputArguments)
%   parsecoastinputs(__, 'Name', value)
%
% Input arguments
%   inputArguments - The input arguments to the function
%       This should usually be varargin itself from the calling function.
%   DefaultUpscale - The fallback value for the upscale argument
%       The default value is 0.
%   DefaultBuffer - The fallback value for the buffer argument
%       The default value is 0.
%   DefaultLatlim - The fallback value for the latlim argument
%       The default value is 90.
%   DefaultMoreBuffers - The fallback value for the morebuffers argument
%       The default value is an empty cell array.
%   DefaultLonOrigin - The fallback value for the lonorigin argument
%       The default value is 180.
%   DefaultForceNew - The fallback value for the forcenew argument
%       The default value is false.
%   DefaultBeQuiet - The fallback value for the bequiet argument
%       The default value is equivalent to 1.
%
% Output arguments
%   upscale - How many times to upscale the data
%   latlim - The inclination angle of the polar caps in degrees
%   buf - The buffer from the coastlines in degrees
%   moreBufs - Additional buffers to apply to the coastlines
%   lonOrigin - The longitude origin of the data
%   forceNew - Force the function to reload the data
%   beQuiet - Suppress the output messages
%       - 0: Show all messages
%       - 1: Suppress all messages not from the function (soft quiet)
%       - 2: Suppress all messages (hard quiet)
%       The default value is 1.
%
% Last modified by 
%   williameclee-at-arizona.edu, 2024/08/07

function varargout = parsecoastinputs(inputArguments, varargin)
    %% Assigning default values
    d = inputParser;
    addRequired(d, 'inputArguments');
    addOptional(d, 'DefaultUpscale', 0);
    addOptional(d, 'DefaultBuffer', 0);
    addOptional(d, 'DefaultLatlim', 90);
    addOptional(d, 'DefaultMoreBuffers', []);
    addOptional(d, 'DefaultLonOrigin', 180);
    addOptional(d, 'DefaultForceNew', false);
    addOptional(d, 'DefaultSaveData', true);
    addOptional(d, 'DefaultBeQuiet', 0.5);
    parse(d, inputArguments, varargin{:});
    upscaleD = d.Results.DefaultUpscale;
    bufD = d.Results.DefaultBuffer;
    latlimD = d.Results.DefaultLatlim;
    moreBufsD = d.Results.DefaultMoreBuffers;
    lonOriginD = d.Results.DefaultLonOrigin;
    forcenewD = d.Results.DefaultForceNew;
    saveDataD = d.Results.DefaultSaveData;
    beQuietD = d.Results.DefaultBeQuiet;

    %% Parsing the 'real' inputs
    p = inputParser;
    p.KeepUnmatched = true;
    addOptional(p, 'Upscale', upscaleD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', bufD, @isnumeric);
    addOptional(p, 'Latlim', latlimD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'MoreBuffers', moreBufsD, @iscell);
    addParameter(p, 'LonOrigin', lonOriginD);
    addParameter(p, 'ForceNew', forcenewD, ...
        @(x) islogical(x) || x == 1 || x == 0);
    addParameter(p, 'SaveData', saveDataD, ...
        @(x) islogical(x));
    addParameter(p, 'BeQuiet', beQuietD, ...
        @(x) islogical(x));
    parse(p, inputArguments{:});

    %% Assigning the parsed values
    upscale = p.Results.Upscale;

    if isempty(upscale) || upscale == 1
        upscale = upscaleD;
    end

    latlim = p.Results.Latlim;

    if any(isempty(latlim)) || any(isnan(latlim))
        latlim = latlimD;
    end

    buf = p.Results.Buffer;

    if isempty(buf)
        buf = bufD;
    end

    moreBufs = p.Results.MoreBuffers;

    if isempty(moreBufs)
        moreBufs = moreBufsD;
    end

    lonOrigin = p.Results.LonOrigin;
    forceNew = logical(p.Results.ForceNew);
    saveData = logical(p.Results.SaveData);
    beQuiet = uint8(p.Results.BeQuiet * 2);

    varargout = ...
        {upscale, latlim, buf, moreBufs, lonOrigin, forceNew, saveData, beQuiet};
end
