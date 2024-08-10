%% OCEANFILENAME
% Return the path to the ocean data file for a given domain.
%
% Syntax
%   dataPath = oceanfilename(domain)
%   dataPath = oceanfilename(domain, 'Name', value)
%   [dataPath, dataFolder, fileExists] = oceanfilename(__)
%
% Input arguments
%   domain - The domain for which to find the ocean data file
%       The domain can be specified as a GeoDomain object, a string, or a
%       cell array of strings.
%       - If a string, the domain is assumed to be the name of the domain.
%       - If a cell array, the first element is the name of the domain, and 
%           the rest are additional arguments.
%   upscale - How many times to upscale the data
%       The default value is 0 (no upscaling).
%   buffer - The (negative) buffer from the coastlines in degrees
%       The default value is 0 (no buffer).
%   latlim - The latitudes of the polar caps in degrees
%       The inclination angle can be specified as a scalar or a 1-by-2
%       vector.
%       - If a scalar, the inclination angle is applied to both the north
%           and south polar caps.
%       - If a 1-by-2 vector, the first element is the inclination angle of
%           the north polar cap, and the second element is the inclination
%           angle of the south polar cap.
%       The default value is 90 (no polar caps).
%   morebuffers - Additional buffers to apply to the coastlines
%       The additional buffers must be specified as a cell array of domain 
%       names and buffer values.
%       The default value is an empty cell array (no additional buffers).
%   RotateBack - A flag to indicate whether the coordinates should be
%       rotated back to their original position.
%       The default value is false (coordinates not rotated).
%
% Output arguments
%   dataPath - The path to the ocean data file
%   dataFolder - The folder containing the ocean data files
%   fileExists - A flag indicating whether the data file exists
%
% Examples
%   dataPath = oceanfilename('indian')
%   dataPath = oceanfilename('indian', 'Upscale', 1, 'Buffer', 1)
%
% Last modified by
%   2024/08/10, williameclee-at-arizona.edu

function [dataPath, dataFolder, fileExists] = ...
        oceanfilename(domain, varargin)
    p = inputParser;
    addRequired(p, 'Domain', ...
        @(x) isa(x, 'GeoDomain') || ischar(x) || iscell(x));
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 0);
    addOptional(p, 'Latlim', 90, ...
        @(x) isnumeric(x) && length(x) <= 2);
    addOptional(p, 'MoreBuffers', []);
    addOptional(p, 'RotateBack', false, ...
        @(x) islogical(x) || isnumeric(x));
    parse(p, domain, varargin{:});
    domain = p.Results.Domain;
    upscale = p.Results.Upscale;
    latlim = p.Results.Latlim;
    buf = p.Results.Buffer;
    moreBuf = p.Results.MoreBuffers;
    rotateBack = p.Results.RotateBack;

    %% Find the data folder
    dataFolder = fullfile(getenv('COASTS'));

    if isempty(dataFolder)
        dataFolder = fullfile(getenv('IFILES'), 'COASTS');
    end

    if ischar(domain)
        domain = GeoDomain(domain, 'Upscale', upscale, 'Buffer', buf, ...
            'Latlim', latlim, 'MoreBuffers', moreBuf);
    elseif iscell(domain)
        domain = GeoDomain(domain{1}, 'Upscale', upscale, ...
            'Buffer', buf, 'Latlim', latlim, 'MoreBuffers', moreBuf, ...
            domain{2:end});
    end

    dataFile = domain.Id;

    if rotateBack
        dataFile = [dataFile, '-rotb'];
    end

    dataPath = fullfile(dataFolder, [dataFile, '.mat']);

    fileExists = exist(dataPath, 'file') == 2;

end
