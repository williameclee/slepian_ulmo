%% AOD1B2PLM
% Auxiliary function to read individual GAC/GAD files
%
% See also
%   AOD1B2PLMT
%
% Authored by
%   2025/03/18, williameclee@arizona.edu (@williameclee)

function varargout = aod1b2plm(varargin)
    %% Initialisation
    % Parse inputs
    [inputPath, Pcenter, smajorAxis, dateFormat] = ...
        parseinputs(varargin);

    % Calculate the midpoint of this data span
    date = finddate(inputPath, dateFormat);

    % Open and scan the file
    [gpSph, gpStdSph] = readaod1bfile(inputPath, smajorAxis);

    % Reorder the coefficients to the ones we use
    if ismember(Pcenter, {'CSR', 'GFZ'})
        gpSph = sortrows(gpSph, [1, 2], 'ascend');
        gpStdSph = sortrows(gpStdSph, [1, 2], 'ascend');
    end

    % TODO: Is WGS84 correction needed?

    varargout = {date, gpSph, gpStdSph};
end

%% Subfunctions
function varargout = parseinputs(inputs)
    dafaultDataCenter = 'CSR';
    defaultSmajorAx = fralmanac('a_EGM96', 'Earth');
    defaultDateFormat = 'datenum';
    ip = inputParser;
    addRequired(ip, 'InputPath', ...
        @(x) exist(x, 'file'));
    addOptional(ip, 'DataCenter', dafaultDataCenter, ...
        @(x) ismember(x, {'CSR', 'GFZ', 'JPL'}) || isempty(x));
    addOptional(ip, 'SemiMajorAxis', defaultSmajorAx, ...
        @(x) isscalar(x) && isnumeric(x) || isempty(x));
    addParameter(ip, 'DateFormat', defaultDateFormat, ...
        @(x) ismember(x, {'datetime', 'datenum'}) || isempty(x));
    parse(ip, inputs{:});
    inputPath = ip.Results.InputPath;
    Pcenter = conddefval(ip.Results.DataCenter, dafaultDataCenter);
    smajorAx = conddefval(ip.Results.SemiMajorAxis, defaultSmajorAx);
    dateFormat = conddefval(ip.Results.DateFormat, defaultDateFormat);

    varargout = {inputPath, Pcenter, smajorAx, dateFormat};
end

function [gpotentialSph, gpotentialStdPlm] = ...
        readaod1bfile(inputPath, smajorAx)
    C = fileread(inputPath);
    C = strsplit(C, '# End of YAML header');
    C = C{2}; % Get the data part
    C = textscan(C, '%s%f%f%f%f%f%f%f%f%s');
    gravitySph = [C{2}, C{3}, C{4}, C{5}];
    gravityStdSph = [C{2}, C{3}, C{6}, C{7}];

    % Convert to geopotential
    gpotentialSph = gravitySph;
    gpotentialSph(:, 3:4) = gpotentialSph(:, 3:4) * smajorAx;
    gpotentialStdPlm = gravityStdSph;
    gpotentialStdPlm(:, 3:4) = gpotentialStdPlm(:, 3:4) * smajorAx;
end

function dateMid = finddate(inputPath, dateFormat)
    [~, fileName] = fileparts(inputPath);
    dateStart = datetime(str2double(fileName(7:10)), ...
        1, str2double(fileName(11:13)));
    dateEnd = datetime(str2double(fileName(15:18)), ...
        1, str2double(fileName(19:21)));
    dateMid = mean([dateStart, dateEnd]);

    if strcmp(dateFormat, 'datenum')
        dateMid = datenum(dateMid); %#ok<DATNM>
    end

end
