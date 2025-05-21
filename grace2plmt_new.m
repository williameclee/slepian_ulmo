%% GRACE2PLMT
% Reads in the Level-2 GRACE geoid products from either the CSR or GFZ data
% centres, does some processing, and saves them as a plmt matrix in a .mat
% file. In particular, the coefficients are reordered to our prefered
% lmcosi format, they are referenced to the WGS84 ellipsoid, the C20/C30
% coefficients are replaced with more accurate measurements from satellite
% laser ranging from Loomis et al, (2020), and the degree one coefficients
% are substituted with those from Sun et al., (2016).  You have the option
% of leaving them as geopotential or converting them to surface mass
% density using the method of Wahr et al. 1998, based on Love numbers.
%
% % Syntax
%   [plmt, dates] = grace2plmt(Pcenter, Rlevel, unit)
%   [plmt, dates] = grace2plmt(__, 'Name', Value)
%
% Input arguments
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data centre is 'CSR'.
%       When the first argument is a cell array, it is interpreted as
%       {Pcenter, Rlevel, Ldata}.
%       Data type: char
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06' (or numbers).
%       The default release level is 'RL06'.
%       Currently, only RL06 is guaranteed to work.
%       Data type: char | [numeric]
%   Ldata - Bandwidth of the date product
%       In the case where there are more than one product from a data
%       centre (such as BA 60 or BB 96 standard L2 products) this allows
%       you to choose between them.
%       The default L is 60.
%       Data type: [numeric]
%   unit - Unit of the output
%       - 'POT': Geopotential field.
%       - 'SD': Surface mass density.
%       The default field is 'SD'.
%       Data type: char
%   TimeRange - Time range of the output
%       When specified, the output will be truncated to the specified time
%       range.
%       The default time range is [] (all available data).
%       Data type: datetime | ([numeric])
%   TimeFormat - Format of the output time
%       - 'datetime': Matlab DATETIME object
%       - 'datenum': Matlab double in DATENUM format
%       The default option is 'datetime'.
%       Data type: char
%   OutputFormat - Format of the output coefficients
%       - 'timefirst': Coefficients are ordered as [time, lmcosi]
%       - 'traditional': Coefficients are ordered as [lmcosi, time]
%       The default option is 'timefirst' (to be consistent with
%       GRACE2PLMT).
%       Data type: char
%   Deg1Correction, C20Correction, C30Correction - Whether to apply these
%       corrections
%       The default options are all true.
%       Data type: logical | ([numeric])
%   ForceNew - Logical flag to force reprocess of the data
%       The default option is false.
%       Data type: logical | ([numeric])
%	SaveData - Logical flag to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%		The default option is true.
%		Data types: logical | ([numeric])
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is false.
%		Data types: logical | ([numeric])
%
% Output arguments
%   plmt - SH coefficients of the gravity field
%       Units: m/s^2 (POT) | kg/m^2 (SD)
%       Data type: DOUBLE
%       Dimension: [nmonths x addmup(L) x 4] | [addmup(L) x 4 x nmonths]
%           (depending on the OutputFormat input)
%   stdPlmt - Standard deviation of the SH coefficients
%       The data type and dimension are the same as for plmt.
%   dates - Time stamps of the coefficients
%       The time stamps are the midpoints of the time intervals of the
%       input data.
%       Datatype: DATATIME | DOUBLE
%           (depending on the TimeFormat input)
%       Dimension: [nmonths x 1]
%
% Data sources
%	GRACE data available from NASA PODAAC:
%	    https://podaac.jpl.nasa.gov/
%
% See also
%   GRACEDEG1, GRACEDEG2, PLM2POT, AOD1B2PLMT
%
% Notes
%   All the intermediate outputs originally printed on the screen are
%   printed to the log file instead.
%
% Last modified by
%   2025/05/21, williameclee@arizona.edu (@williameclee)
%   2022/05/18, charig@email.arizona.edu (@harig00)
%   2020/11/09, lashokkumar@arizona.edu
%   2019/03/18, mlubeck@email.arizona.edu
%   2011/05/17, fjsimons@alum.mit.edu (@fjsimons)

function varargout = grace2plmt_new(varargin)
    %% Initialisation
    % Parse inputs
    [Pcenter, Rlevel, Ldata, unit, timelim, redoDeg1, ...
         deg1corr, c20corr, c30corr, outputFmt, timeFmt, ...
         forceNew, saveData, beQuiet] = ...
        parseinputs(varargin{:});

    % Find the coefficient files
    [inputFolder, outputPath, logPath] = getIOpaths( ...
        Pcenter, Rlevel, Ldata, unit, deg1corr, c20corr, c30corr);

    % If this file already exists, load it.  Otherwise, or if we force it, make
    % a new one (e.g. you added extra months to the database).
    if exist(outputPath, 'file') && ~forceNew
        % Load the peocessed data
        load(outputPath, ...
            'gracePlmt', 'graceStdPlmt', 'dates', 'gravityParam', 'equatorRadius', ...
            'potcoffs', 'date')

        if ~exist('gracePlmt', 'var') && exist('potcoffs', 'var')
            gracePlmt = potcoffs;
        end

        if ~exist('dates', 'var') && exist('date', 'var')
            dates = date;
        end

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), outputPath)
        end

    else
        % Reload from the raw data
        [gracePlmt, graceStdPlmt, dates, gravityParam, equatorRadius] = ...
            grace2plmtCore(Pcenter, Rlevel, Ldata, unit, inputFolder, ...
            deg1corr, c20corr, c30corr, beQuiet, logPath);

        if saveData
            save(outputPath, ...
                'gracePlmt', 'graceStdPlmt', 'dates', 'gravityParam', 'equatorRadius');

            if ~beQuiet
                fprintf('%s saved %s', upper(mfilename), outputPath)
            end

        end

    end

    % Format output
    [gracePlmt, dates] = ...
        formatoutput(gracePlmt, dates, redoDeg1, Pcenter, Rlevel, unit, timelim, outputFmt, timeFmt);

    varargout = {gracePlmt, dates};
end

%% Subfunctions
% Heart of the programme
function [gracePlmt, graceStdPlmt, dates, gravityParam, equatorRadius] = ...
        grace2plmtCore(Pcenter, Rlevel, Ldata, unit, inputFolder, ...
        deg1corr, c20corr, c30corr, beQuiet, logPath)
    %% Computing the coefficients
    logFid = fopen(logPath, 'w');
    fprintf(logFid, '%s - Local time %s\n', ...
        upper(mfilename), datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
    fprintf(logFid, 'Data centre: %s, release level: %s\n', ...
        upper(Pcenter), Rlevel);

    % Find original data files if processed ones are not available
    [inputFiles, Ldata] = getinputfiles(Pcenter, Rlevel, Ldata, inputFolder);
    fprintf(logFid, '%d files to process\n', length(inputFiles));

    % C20 and C30 correction setup
    [tnC2030, tnC2030Std, tnC2030dates] = gracedeg2(Rlevel);

    % Degree 1 correction setup
    [tnDeg1, tnDeg1Std, tnDeg1dates] = gracedeg1(Pcenter, Rlevel);

    % Preallocation
    nDates = length(inputFiles);
    dates = NaT([nDates, 1]);
    gracePlmt = nan([nDates, addmup(Ldata), 4]); % l m cos sin
    graceStdPlmt = nan([nDates, addmup(Ldata), 4]); % l m cos sin

    %% Loop over the months
    wbar = waitbar(0, 'Reading GRACE data', ...
        "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

    for iDate = 1:nDates
        waitbar(iDate / nDates, wbar, ...
            sprintf('Reading GRACE data (%d/%d)', iDate, nDates));

        if getappdata(wbar, 'canceling')
            delete(wbar);
            warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
            'Processing cancelled');
            fclose(logFid);
            return
        end

        % load gravity coefficients
        inputPath = fullfile(inputFolder, inputFiles{iDate});
        [gracePlm, graceStdPlm, date, dateRange] = parsegracesourcefile(inputPath);
        gracePlm(1, 3) = 0;
        dates(iDate) = date;

        fprintf(logFid, '%s - %s\n', date, inputFiles{iDate});

        if deg1corr
            [gracePlm, graceStdPlm] = replaceDeg1(gracePlm, graceStdPlm, date, dateRange, tnDeg1, tnDeg1Std, tnDeg1dates, ...
                beQuiet, logFid);
        end

        if c20corr
            [gracePlm, graceStdPlm] = replaceC20(gracePlm, graceStdPlm, date, dateRange, tnC2030, tnC2030Std, tnC2030dates, ...
                beQuiet, logFid);
        end

        % Now replace C3,0 if possible
        if c30corr
            [gracePlm, graceStdPlm] = replaceC30(gracePlm, graceStdPlm, date, dateRange, tnC2030, tnC2030Std, tnC2030dates, ...
                beQuiet, logFid);
        end

        % Combine into one matrix
        gracePlmt(iDate, :, :) = gracePlm;
        graceStdPlmt(iDate, :, :) = graceStdPlm;
    end

    % WGS84 reference setup
    % For now just hardcode the even zonal coefficients (J), later use
    % Frederik's GRS.m program, don't bother with the higher degrees
    wgs84C20 = 0.108262982131e-2 * -1 / sqrt(5); % will be row 4
    wgs84C40 = -0.237091120053e-5 * -1 / sqrt(5); % will be row 11
    gracePlmt(:, 4, 3) = gracePlmt(:, 4, 3) - wgs84C20;
    gracePlmt(:, 11, 3) = gracePlmt(:, 11, 3) - wgs84C40;

    fclose(logFid);
    delete(wbar);

    %% Converting unit
    % Use the actual parameters stored in the file instead of from FRALMANAC
    [~, ~, ~, ~, gravityParam, equatorRadius] = parsegracesourcefile(inputPath);

    if strcmp(unit, 'POT')
        return
    end

    % Convert gravity to surface geopotential
    gracePlmt(:, :, 3:4) = gracePlmt(:, :, 3:4) * equatorRadius;
    % Convert geopotential to surface mass density
    gracePlmt = plm2pot(gracePlmt, ...
        equatorRadius, gravityParam, equatorRadius, 4);
    graceStdPlmt(:, :, 3:4) = graceStdPlmt(:, :, 3:4) * equatorRadius;
    graceStdPlmt = plm2pot(graceStdPlmt, ...
        equatorRadius, gravityParam, equatorRadius, 4);
end

% Parse input arguments
function varargout = parseinputs(varargin)
    ip = inputParser;
    addOptional(ip, 'Pcenter', 'CSR', ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', 'RL06', ...
        @(x) (isnumeric(x) && isscalar(x)) || ...
        (ischar(validatestring(x, {'RL04', 'RL05', 'RL06'}))));
    addOptional(ip, 'Ldata', 60, ...
        @(x) isnumeric(x) && x > 0);
    addOptional(ip, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'POT', 'SD'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) ((isdatetime(x) || isnumeric(x)) && length(x) == 2) || isempty(x));
    addParameter(ip, 'Deg1Correction', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'C20Correction', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'C30Correction', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'RecomputeDegree1', false, ...
        @(x) islogical(x) || iscell(x) || ischar(x));
    addParameter(ip, 'OutputFormat', 'timefirst', ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'TimeFormat', 'datenum', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));

    if iscell(varargin{1})
        varargin = [varargin{1}{:}, varargin(2:end)];
    end

    parse(ip, varargin{:});

    Pcenter = ip.Results.Pcenter;
    Rlevel = ip.Results.Rlevel;
    Ldata = round(ip.Results.Ldata);
    unit = ip.Results.Unit;
    timelim = ip.Results.TimeRange;
    deg1correction = logical(ip.Results.Deg1Correction);
    c20correction = logical(ip.Results.C20Correction);
    c30correction = logical(ip.Results.C30Correction);
    redoDeg1 = ip.Results.RecomputeDegree1;
    outputFmt = ip.Results.OutputFormat;
    timeFmt = ip.Results.TimeFormat;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);

    if isnumeric(Rlevel)
        Rlevel = sprintf('RL%02d', floor(Rlevel));
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    if islogical(redoDeg1) && redoDeg1
        redoDeg1 = {60, 96, 'ice6gd', GeoDomain('alloceans', "Buffer", 1.5)};
    elseif ischar(redoDeg1)
        redoDeg1 = {60, 96, redoDeg1, GeoDomain('alloceans', "Buffer", 1.5)};
    end

    varargout = ...
        {Pcenter, Rlevel, Ldata, unit, timelim, redoDeg1, ...
         deg1correction, c20correction, c30correction, outputFmt, timeFmt, ...
         forceNew, saveData, beQuiet};
end

% Format the output
function [potcoffs, date] = ...
        formatoutput(potcoffs, date, redoDeg1, Pcenter, Rlevel, unit, timelim, outputFormat, timeFormat)

    if iscell(redoDeg1)
        [~, coeffs, ~] = solvedegree1(Pcenter, Rlevel, redoDeg1{:});
        coeffs = coeffs - mean(coeffs, 1);

        switch unit
            case 'SD' % do nothing
            case 'POT'
                k1 = 0.021; % Swenson et al. 2008
                a = fralmanac('a_EGM96', 'Earth');
                coeffs = coeffs * (1 + k1) / 5517 / a;
        end

        potcoffs(:, 2, 3) = coeffs(:, 1);
        potcoffs(:, 3, 3) = coeffs(:, 2);
        potcoffs(:, 3, 4) = coeffs(:, 3);
    end

    if ~isempty(timelim)
        % Only keep the data within the specified time range
        isValidTime = datetime(date, "ConvertFrom", 'datenum') >= timelim(1) & ...
            datetime(date, "ConvertFrom", 'datenum') <= timelim(2);
        potcoffs = potcoffs(isValidTime, :, :);
        date = date(isValidTime);
    end

    switch outputFormat
        case 'timefirst'
            % Do nothing
        case 'traditional'
            potcoffs = permute(potcoffs, [2, 3, 1]);
    end

    switch timeFormat
        case 'datenum'

            if isdatetime(date)
                date = datenum(date); %#ok<DATNM>
            end

        case 'datetime'

            if isnumeric(date)
                date = datetime(date, "ConvertFrom", 'datenum');
            end

    end

end

% Get the input folder and output file names
function [inputFolder, outputPath, logPath] = ...
        getIOpaths(Pcenter, Rlevel, Ldata, unit, deg1corr, c20corr, c30corr)

    if ~isempty(getenv('ORIGINALGRACEDATA'))
        inputFolder = fullfile(getenv('ORIGINALGRACEDATA'), ...
            Rlevel, Pcenter);
    elseif ~isempty(getenv('GRACEDATA'))
        inputFolder = fullfile(getenv('GRACEDATA'), 'raw', ...
            Rlevel, Pcenter);
    else
        inputFolder = fullfile(getenv('IFILES'), 'GRACE', 'raw', ...
            Rlevel, Pcenter);
    end

    % Where you would like to save the new .mat file
    if ~isempty(getenv('GRACEDATA'))
        outputFolder = fullfile(getenv('GRACEDATA'));
    else
        outputFolder = fullfile(getenv('IFILES'), 'GRACE');
    end

    switch unit % no otherwise case since input validity is already checked
        case 'SD'
            outputFile = sprintf('%s_%s_alldata_%s_%s.mat', ...
                Pcenter, Rlevel, num2str(Ldata), unit);
        case 'POT'
            outputFile = sprintf('%s_%s_alldata_%s.mat', ...
                Pcenter, Rlevel, num2str(Ldata));
    end

    if ~c30corr
        outputFile = strrep(outputFile, 'alldata', 'alldata_nC30');
    end

    if ~c20corr
        outputFile = strrep(outputFile, 'alldata', 'alldata_nC20');
    end

    if ~deg1corr
        outputFile = strrep(outputFile, 'alldata', 'alldata_nDeg1');
    end

    outputPath = fullfile(outputFolder, outputFile);
    logPath = fullfile(outputFolder, sprintf('%s_log.txt', upper(mfilename)));

    if ~isfile(outputPath) == 2 && ~exist(inputFolder, 'dir') == 2
        error('The data you asked for are not currently stored\nPlease check the input folder %s', inputFolder)
    end

end

% Get the raw input files
function [dataFiles, Ldata] = ...
        getinputfiles(Pcenter, Rlevel, Ldata, inputFolder)
    % Only RL06 is supported for now
    if ~strcmp(Rlevel, 'RL06')
        error( ...
            sprintf('%s:LoadData:SolutionLoadingNotImplemented', upper(mfilename)), ...
            'Loading %s %s solutions is not currently implemented', ...
            upper(Pcenter), Rlevel);
    end

    % Get the data files
    switch Ldata
        case 60
            dataFiles = ls2cell(fullfile(inputFolder, ...
            'GSM-2_*_BA01_06*'));
        case 96
            dataFiles = ls2cell(fullfile(inputFolder, ...
            'GSM-2_*_BB01_06*'));
        otherwise
            error(sprintf('%s:LoadData:SolutionDNE', upper(mfilename)), ...
                '%s degree %d solutions not available', ...
                upper(Pcenter), Ldata);
    end

    % Make sure the files exist
    if isempty(dataFiles)
        error(sprintf('%s:LoadData:NoRawGRACEDataFound', upper(mfilename)), ...
            'No data files found in %s', inputFolder)
    end

end

% Do degree 1 correction
function [gracePlm, graceStdPlm] = ...
        replaceDeg1(gracePlm, graceStdPlm, date, dateRange, tnDeg1, tnDeg1Std, tnDeg1dates, beQuiet, fid)
    % Find a degree 1 data point that is within the month of our GRACE data
    iTn = tnDeg1dates(:) > dateRange(1) & tnDeg1dates(:) < dateRange(2) & ...
        ~any(isnan(tnDeg1(:, :, 3:4)), [2, 3]);

    if ~any(iTn)

        fprintf(fid, '  C10: (not replaced)   %+10.6e\n', ...
            gracePlm(2, 3));

        if ~beQuiet
            warning(sprintf('%s:NoTNReplacementAvailable', upper(mfilename)), ...
                'No degree 1 replacement for %s', date);
        end

        return

    end

    if sum(iTn) > 1
        % We have more than one month. Use the closest value
        [~, iTn] = min(abs(date - tnDeg1dates));
    end

    fprintf(fid, '  C10: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(2, 3), tnDeg1(iTn, 1, 3));
    fprintf(fid, '  C11: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(3, 3), tnDeg1(iTn, 2, 3));
    fprintf(fid, '  S11: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(3, 4), tnDeg1(iTn, 2, 4));

    % Replacement
    gracePlm(2:3, 3:4) = squeeze(tnDeg1(iTn, :, 3:4));
    graceStdPlm(2:3, 3:4) = squeeze(tnDeg1Std(iTn, :, 3:4));
end

% Do C20 correction
function [gracePlm, graceStdPlm] = ...
        replaceC20(gracePlm, graceStdPlm, date, dateRange, tnC2030, tnC2030Std, tnC2030dates, beQuiet, fid)
    iTn = tnC2030dates > dateRange(1) & tnC2030dates < dateRange(2) & ...
        ~isnan(tnC2030(:, 1));

    if ~any(iTn)

        fprintf(fid, '  C20: (not replaced)   %+10.6e\n', ...
            gracePlm(4, 3));

        if ~beQuiet
            warning(sprintf('%s:NoTNReplacementAvailable', upper(mfilename)), ...
                'No C20 replacement for %s', date);
        end

        return
    end

    if sum(iTn) > 1
        % We have more than one month. Use the closest value
        [~, iTn] = min(abs(date - tnC2030dates));
    end

    fprintf(fid, '  C20: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(4, 3), tnC2030(iTn, 1));

    % Replacement
    gracePlm(4, 3) = tnC2030(iTn, 1);
    graceStdPlm(4, 3) = tnC2030Std(iTn, 1);
end

% Do C20 correction
function [gracePlm, graceStdPlm] = ...
        replaceC30(gracePlm, graceStdPlm, date, dateRange, tnC2030, tnC2030Std, tnC2030dates, beQuiet, fid)
    iTn = tnC2030dates > dateRange(1) & tnC2030dates < dateRange(2) & ...
        ~isnan(tnC2030(:, 2));

    if ~any(iTn)

        fprintf(fid, '  C30: (not replaced)   %+10.6e\n', ...
            gracePlm(7, 3));

        if ~beQuiet
            warning(sprintf('%s:NoTNReplacementAvailable', upper(mfilename)), ...
                'No C30 replacement for %s', date);
        end

        return
    end

    if sum(iTn) > 1
        % We have more than one month. Use the closest value
        [~, iTn] = min(abs(date - tnC2030dates));
    end

    fprintf(fid, '  C30: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(7, 3), tnC2030(iTn, 2));

    % Replacement
    gracePlm(7, 3) = tnC2030(iTn, 2);
    graceStdPlm(7, 3) = tnC2030Std(iTn, 2);
end
