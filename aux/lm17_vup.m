%% LM17_VUP
% Processes the LM17.3 GIA model to spherical harmonics vertical uplift.
%
% Last modified by
%   2024/10/03, williameclee@arizona.edu (@williameclee)

function lmcosiM = lm17_vup(varargin)
    %% Initialisation
    % Parse inputs
    ip = inputParser;
    addOptional(ip, 'inputFolder', '', @ischar);
    addOptional(ip, 'outputPath', '', @ischar);
    parse(ip, varargin{:});
    inputFolder = ip.Results.inputFolder;
    outputPath = ip.Results.outputPath;

    % Locate I/O
    if isempty(inputFolder)
        inputFolder = fullfile(getenv('IFILES'), 'GIA', 'LM17.3');
    end

    if isempty(outputPath)
        outputPath = fullfile(inputFolder, 'LM17.3_vup.mat');
    end

    if exist(outputPath, 'file')
        load(outputPath, 'lmcosiM')

        if exist('lmcosiM', 'var')
            fprintf('%s loaded %s\n', upper(mfilename), outputPath)
            return
        end

    end

    inputPath = fullfile(inputFolder, 'LM17.3_0.5x0.5_uplift_globe.txt');

    if ~isfile(inputPath)
        error('File not found: %s\nPlease download it from %s', ...
            inputPath, 'https://sites.google.com/view/holgersteffenlm/startseite/data')
    end

    %% Load data
    T = readtable(inputPath);
    T.Long = mod(T.Long, 360);

    %% Process data
    % Interpolation
    h = 0.5;
    lon = 0:h:360;
    lat = flip(-90:h:90);
    [lonn, latt] = meshgrid(lon, lat);
    F = scatteredInterpolant(T.Long, T.Lat, T{:, 3}, 'linear', 'nearest');
    z = F(lonn, latt);

    % Convert to spherical harmonics surface density
    wPlm = xyz2plm(z, 180);

    %% Save data
    lmcosiM = wPlm;
    save(outputPath, 'lmcosiM')
    fprintf('%s saved %s\n', upper(mfilename), outputPath)
end
