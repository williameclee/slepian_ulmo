function varargout = closecoastline(varargin)

    if nargin == 1
        lonlat = varargin{1};
    elseif nargin == 2
        lonlat = [varargin{1}, varargin{2}];
    end

    % remove coordinates with ONLY ONE NaN
    lonlat = lonlat(all(isnan(lonlat), 2) | all(~isnan(lonlat), 2), :);

    % Split coastlines into separate segments based on NaNs
    lonlatSplit = splitxy(lonlat);

    % Close each segment
    for iCell = 1:length(lonlatSplit)

        if any(lonlatSplit{iCell}(1, :) ~= lonlatSplit{iCell}(end, :))
            lonlatSplit{iCell} = ...
                [lonlatSplit{iCell}; lonlatSplit{iCell}(1, :)];
        end

    end

    % Join the segments back together
    lonlatClosed = joinxy(lonlatSplit);

    if nargout == 1
        varargout{1} = lonlatClosed;
    elseif nargout == 2
        varargout{1} = lonlatClosed(:, 1);
        varargout{2} = lonlatClosed(:, 2);
    end

end
