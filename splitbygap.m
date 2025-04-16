%% SPLITBYGAP
% Inserts NaN (or NaT) into a time series to split it at gaps greater than a given threshold.
%
% Syntax
% [time, series1, series2, ...] = splitbygap(time, tol, series1, series2, ...)
%
% Inputs
%   time - A DATETIME, DURATION, or numeric vector
%   tol - A scalar tolerance value
%       If empty, defaults to 120 days for DATETIME/DURATION or 1 for numeric.
%   series1, series2, ... - Vectors or matrices of the same length as time
%       These will be split at the same indices as time.
%
% Outputs
%   time - The modified time vector with NaN (or NaT) inserted
%   series1, series2, ... - The modified series vectors/matrices with NaN inserted
%
% Last modified by
%   2025/04/16, williameclee@arizona.edu (@williameclee)

function [x, varargout] = splitbygap(x, tol, varargin)

    if isempty(tol)

        if isdatetime(x) || isduration(x)
            tol = days(120);
        else
            tol = 1;
        end

    end

    if size(x, 1) == 1
        x = x(:);
        transX = true;
    else
        transX = false;
    end

    transY = cell([length(varargin), 1]);
    padY = ones([length(varargin), 1]);

    for i = 1:length(varargin)

        if size(varargin{i}, 1) == 1
            varargin{i} = varargin{i}';
            transY{i} = true;
        else
            transY{i} = false;
        end

        padY(i) = size(varargin{i}, 2);

    end

    if isdatetime(x) || isduration(x)
        splitObj = NaT;
    else
        splitObj = nan;
    end

    xGap = diff(x);
    % Split the time series if there is a gap
    splitIndex = find(xGap > tol);

    if ~isempty(splitIndex)

        for iSplit = length(splitIndex):-1:1
            x = ...
                [x(1:splitIndex(iSplit)); splitObj; ...
                 x(splitIndex(iSplit) + 1:end)];

            for i = 1:length(varargin)
                varargin{i} = ...
                    [varargin{i}(1:splitIndex(iSplit), :); ...
                     nan([1, padY(i)]); ...
                     varargin{i}(splitIndex(iSplit) + 1:end, :)];
            end

        end

    end

    if transX
        x = x';
    end

    for i = 1:length(varargin)

        if transY{i}
            varargin{i} = varargin{i}';
        end

    end

    varargout = varargin;

end
