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
    parse(p, domain, varargin{:});
    domain = p.Results.Domain;
    upscale = p.Results.Upscale;
    latlim = p.Results.Latlim;
    buf = p.Results.Buffer;
    moreBuf = p.Results.MoreBuffers;

    %% Find the data folder
    dataFolder = fullfile(getenv('COASTS'));

    if isempty(dataFolder)
        dataFolder = fullfile(getenv('IFILES'), 'COASTS');
    end

    if isa(domain, 'GeoDomain')
        dataFile = [domain.Id, '.mat'];
    else
        dataFileAttr = dataattrchar('Upscale', upscale, 'Buffer', buf, ...
            'Latlim', latlim, 'MoreBuffers', moreBuf);

        dataFile = [capitalise(domain), '-', dataFileAttr, '.mat'];
    end

    dataPath = fullfile(dataFolder, dataFile);

    fileExists = exist(dataPath, 'file') == 2;

end
