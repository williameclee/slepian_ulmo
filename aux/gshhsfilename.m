function [file, folder, fileExists] = gshhsfilename(varargin)
    %% Initialisation
    p = inputParser;
    addOptional(p, 'DataQuality', 'c', ...
        @(x) ischar(x) && ismember(x, 'cfhil'));
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'MinLandArea', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Tolerence', 0, ...
        @(x) isnumeric(x) || isempty(x));
    parse(p, varargin{:});

    dataQuality = p.Results.DataQuality;
    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    minLandArea = p.Results.MinLandArea;
    tol = p.Results.Tolerence;

    %% Generating file name
    if upscale == 0 || upscale == 1 || isempty(upscale)
        upscaleString = '';
    else
        upscaleString = ['-', num2str(upscale)];
    end

    if buf == 0 || isempty(buf)
        bufString = '';
    else
        bufString = ['-', num2str(buf)];

        if isempty(upscaleString)
            upscaleString = '-0';
        end

    end

    if minLandArea == 0 || isempty(minLandArea)
        minLandAreaString = '';
    else
        minLandAreaString = ['-', num2str(minLandArea)];

        if isempty(upscaleString)
            upscaleString = '-0';
        end

        if isempty(bufString)
            bufString = '-0';
        end

    end

    if tol == 0 || isempty(tol)
        tolString = '';
    else
        tolString = ['-', num2str(tol)];

        if isempty(upscaleString)
            upscaleString = '-0';
        end

        if isempty(bufString)
            bufString = '-0';
        end

        if isempty(minLandAreaString)
            minLandAreaString = '-0';
        end

    end

    folder = getenv('GSHHS');
    file = fullfile(folder, ...
        ['gshhs_', dataQuality, upscaleString, bufString, ...
         minLandAreaString, tolString, '.mat']);
    fileExists = exist(file, 'file') == 2;
end
