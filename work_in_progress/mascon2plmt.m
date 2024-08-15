function varargout = mascon2plmt(varargin)
    p = inputParser;
    addOptional(p, 'dataType', 'mascon', @ischar);
    addOptional(p, 'L', 60, @isnumeric);
    addOptional(p, 'meshSize', 0.5, @isnumeric);
    parse(p, varargin{:});

    dataType = p.Results.dataType;
    L = p.Results.L;
    meshSize = p.Results.meshSize;

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

    if exist(outputPath, 'file')
        load(outputPath, 'time', 'dataPlmt');

        varargout = {dataPlmt, time};
        return;
    end

    time = ncread(dataPath, 'time');
    time = datetime(2002, 1, 1) + double(time) - 1;
    lat = ncread(dataPath, 'lat');
    lon = ncread(dataPath, 'lon');
    dataF = ncread(dataPath, 'lwe_thickness');

    dataPlmt = zeros([length(time), addmup(L), 4]);

    parfor t = 1:length(time)
        dataPlmt(t, :, :) = mascon2plm(dataF(:, :, t), lon, lat, L, meshSize);
    end

    dataPlmt(:, :, 3:4) = dataPlmt(:, :, 3:4) * 10; % From cm to kg/m^2

    fprintf('%s saving %s\n', upper(mfilename), outputPath);
    save(outputPath, 'time', 'dataPlmt');

    varargout = {dataPlmt, time};
end

function masconPlm = mascon2plm(mascon, lon, lat, L, meshSize)
    [lonn, latt] = ndgrid(lon, lat);
    newLon = 0:meshSize:360;
    newLat = -90:meshSize:90;
    [newLonn, newLatt] = ndgrid(newLon, newLat);
    F = griddedInterpolant(lonn, latt, double(mascon), 'linear', 'linear');
    masconC = F(newLonn, newLatt)';
    masconPlm = xyz2plm(flip(masconC), L);
end
