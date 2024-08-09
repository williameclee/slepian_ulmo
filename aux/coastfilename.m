function [dataPath, dataFolder, fileExists] = ...
        coastfilename(domain, varargin)
    p = inputParser;
    addRequired(p, 'Domain');
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

    if isa(domain, 'GeoDomain')
        dataFile = domain.Id;
    else
        dataFileAttr = dataattrchar('Upscale', upscale, 'Buffer', buf, ...
            'Latlim', latlim, 'MoreBuffers', moreBuf);

        dataFile = [capitalise(domain), dataFileAttr];
    end

    if rotateBack
        dataFile = [dataFile, '-1'];
    end

    dataPath = fullfile(dataFolder, [dataFile, '.mat']);

    fileExists = exist(dataPath, 'file') == 2;

end
