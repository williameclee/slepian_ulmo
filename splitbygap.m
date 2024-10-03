function [x, varargout] = splitbygap(x, tol, varargin)

    if isempty(tol)

        if isdatetime(x)
            tol = days(120);
			splitObj = NaT;
        else
            tol = 1;
			splitObj = nan;
        end

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
