%% PLOTDQM
% Plots a quick and dirty map (on a normal axes) of the given data.
%
% Syntax
%   plotqdm(lonlat)
%   plotqdm(lon, lat)
%   plotqdm(__, c)
%   plotqdm(__, 'Name', Value)
%   h = plotqdm(__)
%
% Input arguments
%   lonlat - 2-column matrix of longitude and latitude.
%   lon, lat - Vectors of longitude and latitude.
%   c - Colour of the line.
%   The name-value pair arguments are passed to the plot function.
%
% Output arguments
%   h - Handle to the axes.
%
% Last modified by
%   2024/08/12, williameclee@arizona.edu (@williameclee)

function h = plotqdm(varargin)
    [lon, lat, colour, Inputs] = parseinputs(varargin);

    lon = lon - floor(min(lon) / 360) * 360;
    xLim = [min(lon), max(lon)] + [-5, 5];

    if xLim(2) - xLim(1) == 370
        xLim = [0, 360];
    end

    yLim = [min(lat), max(lat)] + [-5, 5];
    yLim = min(max(yLim, -90), 90);

    %% Plotting
    plot(lon, lat, 'Color', colour, "HandleVisibility", 'off', ...
        Inputs{:})
    axis equal
    axis tight
    grid on

    xlim(xLim)
    ylim(yLim)

    xticklabels(cellfun(@formatlonticks, num2cell(xticks), ...
        "UniformOutput", false))
    yticklabels(cellfun(@formatlatticks, num2cell(yticks), ...
        "UniformOutput", false))

    set(gca, 'Box', 'on', 'Layer', 'top')

    if nargout == 0
        return
    end

    h = gca;
end

%% Subfunctions
function [lon, lat, c, Inputs] = parseinputs(Inputs)
    startOfOptional = find(cellfun(@ischar, Inputs), 1);

    if ~isempty(startOfOptional) && startOfOptional == 1
        error('The first argument(s) must be the data.');
    end

    if startOfOptional >= 5
        error('Invalid data arguments.');
    elseif isempty(startOfOptional)
        startOfOptional = length(Inputs) + 1;
    end

    if startOfOptional == 2

        if ~any(size(Inputs{1}) == 2)
            error('The first argument must be a 2-column matrix.');
        end

        if size(Inputs{1}, 2) == 2
            lon = Inputs{1}(:, 1);
            lat = Inputs{1}(:, 2);
        else
            lon = Inputs{1}(1, :)';
            lat = Inputs{1}(2, :)';
        end

        Inputs = Inputs(2:end);

    elseif startOfOptional >= 3
        lon = Inputs{1};
        lat = Inputs{2};
        Inputs = Inputs(3:end);
    end

    lon = lon(:);
    lat = lat(:);

    if ~isequal(size(lon), size(lat))
        error('Longitude and latitude must have the same size.');
    end

    c = 'k';

    if isempty(Inputs)
        return
    end

    if ~(ischar(Inputs{1}) || isstring(Inputs{1})) && ...
            ~(isnumeric(Inputs{1}) && length(Inputs{1}) == 3)
        error('Unrecognised colour.');
    elseif ischar(Inputs{1}) || isstring(Inputs{1})

        if ismember(Inputs{1}, ...
                ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w', ...
                 'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black', 'white'])
            c = Inputs{1};
            Inputs = Inputs(2:end);
        end

    else
        c = Inputs{1};
        Inputs = Inputs(2:end);
    end

end
