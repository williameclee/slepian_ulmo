%% PARSEGRACESOURCEFILE
% Parses a GRACE source file and returns the gravity field and its uncertainty. Should work for both GSM and GAC/GAD files.
%
% See also
%	GRACE2PLMT (GRACE2PLMT_NEW), AOD1B2PLMT
%
% Notes
%	This is a helper function with limited documentation.
%
% Authored by
%   2025/05/20, williameclee@arizona.edu (@williameclee)
function varargout = parsegracesourcefile(dataPath)
    %% Loading data
    % Extract the data from the file
    dataBarrier = '# End of YAML header\n';
    data = fileread(dataPath);
    data = strsplit(data, dataBarrier);
    data = data{2};
    data = textscan(data, '%s%f%f%f%f%f%f%s%s%s');
    gravitySph = [data{2}, data{3}, data{4}, data{5}];
    gravityStdSph = [data{2}, data{3}, data{6}, data{7}];

    % Sort the data
    gravitySph = sortrows(gravitySph, [1, 2], 'ascend');
    gravityStdSph = sortrows(gravityStdSph, [1, 2], 'ascend');

    % Pad the degree 0 and 1 data (for JPL)
    if gravitySph(1, 1) == 2
        gravitySph = ...
            [[0, 0, 1, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravitySph];
    end

    if gravityStdSph(1, 1) == 2
        gravityStdSph = ...
            [[0, 0, 1, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravityStdSph];
    end

    %% Quality check
    % Make sure the sizes of the matrices are correct
    if ~isequal(size(gravitySph), size(gravityStdSph))
        error(sprintf('%s:DataSizeNotMatch', upper(mfilename)), ...
            'Size of the gravity field (%d) and the uncertainty (%d) does not match. Please check the input file %s', ...
            size(gravitySph, 1), size(gravityStdSph, 1), dataPath)
    end

    if size(gravitySph, 1) ~= addmup(gravitySph(end, 1))
        error(sprintf('%s:DataSizeNotMatch', upper(mfilename)), ...
            'Size of the gravity field (%d) does not match the expectation (%d) for degree %d', ...
            size(gravitySph, 1), addmup(gravitySph(end, 1)), gravitySph(end, 1))
    end

    %% Finalising output
    if nargout == 2
        varargout = {gravitySph, gravityStdSph};
        return
    end

    startDate = datetime(data{8}{1}, ...
        "InputFormat", 'yyyyMMdd.hhmm', "Format", 'yyyy/MM/dd HH:mm');
    endDate = datetime(data{9}{1}, ...
        "InputFormat", 'yyyyMMdd.hhmm', "Format", 'yyyy/MM/dd HH:mm');
    meanDate = mean([startDate, endDate]);

    varargout = {gravitySph, gravityStdSph, meanDate, [startDate, endDate]};
end
