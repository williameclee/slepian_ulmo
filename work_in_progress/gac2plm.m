function varargout = gac2plm(varargin)
    %% Initialisation
    % Parse inputs
    [inputPath, Pcenter, Ldata, a, unit, dateFormat] = ...
        parseinputs(varargin);

    % Calculate the midpoint of this data span
    date = finddate(inputPath, dateFormat);

    % Open and scan the file
    [gpotentialPlm, gpotentialStdPlm] = readgacfile(inputPath, a);

    if ismember(Pcenter, {'CSR', 'GFZ'})
        gpotentialPlm = reorderplmcoeffs(gpotentialPlm, Ldata);
        gpotentialStdPlm = reorderplmcoeffs(gpotentialStdPlm, Ldata);
    end

    % TODO: Is WGS84 correction needed?

    if strcmp(unit, 'POT')
        varargout = {date, gpotentialPlm, gpotentialStdPlm};
        return
    end

    % Convert the geopotential coefficients into surface mass density,
    sdensityPlm = plm2pot(gpotentialPlm, [], [], [], 4);

    if nargout == 2
        varargout = {date, sdensityPlm};
        return
    end

    sdensityStdPlm = plm2pot(gpotentialStdPlm, [], [], [], 4);

    varargout = {date, sdensityPlm, sdensityStdPlm};
end

%% Subfunctions
function varargout = parseinputs(inputs)
    dafaultDataCenter = 'CSR';
    defaultLdata = 180;
    defaultSmajorAx = fralmanac('a_EGM96', 'Earth');
    defaultUnit = 'SD';
    defaultDateFormat = 'datenum';
    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p, 'InputPath', ...
        @(x) exist(x, 'file'));
    addOptional(p, 'DataCenter', dafaultDataCenter, ...
        @(x) ismember(x, {'CSR', 'GFZ', 'JPL'}) || isempty(x));
    addOptional(p, 'Ldata', defaultLdata, ...
        @(x) isscalar(x) && isnumeric(x));
    addOptional(p, 'SemiMajorAxis', defaultSmajorAx, ...
        @(x) isscalar(x) && isnumeric(x) || isempty(x));
    addOptional(p, 'Unit', defaultUnit, ...
        @(x) ismember(x, {'POT', 'SD'}) || isempty(x));
    addParameter(p, 'DateFormat', defaultDateFormat, ...
        @(x) ismember(x, {'datetime', 'datenum'}) || isempty(x));
    parse(p, inputs{:});
    inputPath = p.Results.InputPath;
    Pcenter = conddefval(p.Results.DataCenter, dafaultDataCenter);
    Ldata = conddefval(p.Results.Ldata, defaultLdata);
    smajorAx = conddefval(p.Results.SemiMajorAxis, defaultSmajorAx);
    unit = conddefval(p.Results.Unit, defaultUnit);
    dateFormat = conddefval(p.Results.DateFormat, defaultDateFormat);

    varargout = {inputPath, Pcenter, Ldata, smajorAx, unit, dateFormat};
end

function [gpotentialPlm, gpotentialStdPlm] = ...
        readgacfile(inputPath, smajorAx)
    fid = fopen(inputPath);
    C = textscan(fid, '%s%s%s%s%s%s%s%s%s%s');
    fclose(fid);

    isData = cellfun(@(x) strcmp(x, 'GRCOF2'), C{1});
    gravityPlm = [C{2}(isData), C{3}(isData), C{4}(isData), C{5}(isData)];
    gravityPlm = cellfun(@str2double, gravityPlm);
    gravityStdPlm = [C{2}(isData), C{3}(isData), C{6}(isData), C{7}(isData)];
    gravityStdPlm = cellfun(@str2double, gravityStdPlm);

    % Convert to geopotential
    gpotentialPlm = gravityPlm;
    gpotentialPlm(:, 3:4) = gpotentialPlm(:, 3:4) * smajorAx;
    gpotentialStdPlm = gravityStdPlm;
    gpotentialStdPlm(:, 3:4) = gpotentialStdPlm(:, 3:4) * smajorAx;
end

function plmCoeffs_reordered = reorderplmcoeffs(plmCoeffs, Ldata)
    [dems, dels] = addmon(Ldata);

    plmCoeffs_reordered = zeros(size(plmCoeffs));
    revdel = [0, Ldata:-1:0];
    i = 1;

    for j = 1:length(dems)
        k = dels(i) + 1 + sum(revdel((1:dems(i) + 1)));
        plmCoeffs_reordered(j, :) = plmCoeffs(k, :);
        i = i + 1;
    end

end

function dateMid = finddate(inputPath, dateFormat)
    %#ok<*DATNM>
    [~, fileName] = fileparts(inputPath);
    dateStart = datetime(str2double(fileName(7:10)), ...
        1, str2double(fileName(11:13)));
    dateStart = datenum(dateStart);
    dateEnd = datetime(str2double(fileName(15:18)), ...
        1, str2double(fileName(19:21)));
    dateEnd = datenum(dateEnd);
    dateMid = (dateStart + dateEnd) / 2;

    if strcmp(dateFormat, 'datetime')
        dateMid = datetime(dateMid, 'ConvertFrom', 'datenum');
    end

end
