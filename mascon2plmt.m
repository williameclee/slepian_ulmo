%% MASCON2PLMT
% Convert (mesh) mascon data to spherical harmonics
%
% Last modified by
%   2024/08/30, williameclee@arizona.edu (@williameclee)

function varargout = mascon2plmt(varargin)
    p = inputParser;
    addOptional(p, 'dataType', 'mascon', @ischar);
    addOptional(p, 'L', 60, @isnumeric);
    addOptional(p, 'meshSize', 0.5, @isnumeric);
    addParameter(p, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    dataType = p.Results.dataType;
    L = p.Results.L;
    meshSize = p.Results.meshSize;
    forceNew = p.Results.ForceNew;
    beQuiet = p.Results.BeQuiet;

    dataFolder = fullfile(getenv('GRACEDATA'), 'mascon', 'CSR', 'RL06');

    switch dataType
        case 'mascon'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_all-corrections.nc';
        case 'gsu'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_GSU-component.nc';
        case 'gia'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_GIA-component.nc';
        case 'gad'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_GAD-component.nc';
        case 'masconC20'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_MasconC20-component.nc';
        case 'masconC30'
            dataFile = 'CSR_GRACE-FO_RL0602_Mascons_MasconC30-component.nc';
        case 'slrC20'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_SLR-C20-component.nc';
        case 'slrC30'
            dataFile = 'CSR_GRACE-FO_RL0602_Mascons_SLR-C30-component.nc';
        case 'deg1'
            dataFile = 'CSR_GRACE_GRACE-FO_RL0602_Mascons_degree1-component.nc';
        otherwise
            error('Invalid data type');
    end

    dataPath = fullfile(dataFolder, dataFile);

    switch dataType
        case 'mascon'
            outputFile = sprintf('CSR-RL06-%i-SD-mascon.mat', L);
        otherwise
            outputFile = sprintf('CSR-RL06-%i-SD-mascon-%s.mat', L, dataType);
    end

    outputPath = fullfile(dataFolder, outputFile);

    if exist(outputPath, 'file') && ~forceNew
        load(outputPath, 'time', 'dataPlmt');

        if size(dataPlmt, 2) ~= 4
            dataPlmt = permute(dataPlmt, [2, 3, 1]);
            save(outputPath, 'time', 'dataPlmt');

            if ~beQuiet
                fprintf('%s loaded and updated %s\n', upper(mfilename), outputPath);
            end

        else

            if ~beQuiet
                fprintf('%s loaded %s\n', upper(mfilename), outputPath);
            end

        end

        varargout = {dataPlmt, time};
        return
    end

    time = ncread(dataPath, 'time');
    time = datetime(2002, 1, 0) + double(time);
    lat = ncread(dataPath, 'lat');
    lon = ncread(dataPath, 'lon');
    dataF = ncread(dataPath, 'lwe_thickness');

    dataPlmt = zeros([addmup(L), 4, length(time)]);

    parfor t = 1:length(time)
        dataPlmt(:, :, t) = ...
            mascon2plm(dataF(:, :, t), lon, lat, L, meshSize);
    end

    dataPlmt(:, 3:4, :) = dataPlmt(:, 3:4, :) * 10; % From cm to kg/m^2

    save(outputPath, 'time', 'dataPlmt');

    if ~beQuiet
        fprintf('%s saved %s\n', upper(mfilename), outputPath);
    end

    varargout = {dataPlmt, time};
end

function masconPlm = mascon2plm(mascon, lon, lat, L, ~)
    [lonn, latt] = ndgrid(lon, lat);
    meshSize = abs(lon(2) - lon(1));
    newLon = 0:meshSize:360;
    newLat = -90:meshSize:90;
    [newLonn, newLatt] = ndgrid(newLon, newLat);
    F = griddedInterpolant(lonn, latt, double(mascon), 'linear', 'linear');
    masconC = F(newLonn, newLatt)';
    masconPlm = xyz2plm(flip(masconC'), L);
end
