%% SPLITBYGAP
% Inserts NaN (or NaT) into a time series to split it at gaps greater than a given threshold.
%
% Syntax
% [time, series1, series2, ...] = splitbygap(time, tol, series1, series2, ...)
%
% Last modified by
%   2025/03/28, williameclee@arizona.edu (@williameclee)

function [x, varargout] = splitbygap(x, tol, varargin)

    if isempty(tol)

        if isdatetime(x) || isduration(x)
            tol = days(120);
        else
            tol = 1;
        end

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
                [x(1:splitIndex(iSplit)), splitObj, ...
                 x(splitIndex(iSplit) + 1:end)];

            for i = 1:length(varargin)
                varargin{i} = ...
                    [varargin{i}(1:splitIndex(iSplit)), nan, ...
                     varargin{i}(splitIndex(iSplit) + 1:end)];
            end

        end

    end

    varargout = varargin;

end
