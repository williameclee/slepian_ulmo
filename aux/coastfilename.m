function [dataFile, dataFolder, fileExists] = ...
        coastfilename(region, varargin)
    p = inputParser;
    addRequired(p, 'region');
    addOptional(p, 'upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'buffer', 0);
    addOptional(p, 'inclang', 90, ...
        @(x) isnumeric(x) && length(x) <= 2);
    addOptional(p, 'MoreBuffer', []);
    parse(p, region, varargin{:});
    region = p.Results.region;
    upscale = p.Results.upscale;
    inclang = p.Results.inclang;
    buf = p.Results.buffer;
    moreBuf = p.Results.MoreBuffer;

    %% Find the data folder
    dataFolder = fullfile(getenv('COASTS'));

    if isempty(dataFolder)
        dataFolder = fullfile(getenv('IFILES'), 'COASTS');
    end

    dataFileAttr = dataattrchar('Upscale', upscale, 'Buffer', buf, 'Inclang', inclang, 'MoreBuffer', moreBuf);

    dataFile = [capitalise(region), '-', dataFileAttr, '.mat'];
    dataFile = fullfile(dataFolder, dataFile);

    fileExists = exist(dataFile, 'file') == 2;

end
