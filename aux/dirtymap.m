function dirtymap(varargin)

    if nargin == 1
        lon = varargin{1}(:, 1);
        lat = varargin{1}(:, 2);
    elseif nargin == 2
        lon = varargin{1};
        lat = varargin{2};
    else
        error('Invalid number of arguments');
    end

    plot(lon, lat, 'k', "HandleVisibility", "off")
    axis equal
    axis tight
    grid on
    xlim([min(lon), max(lon)] + [-5, 5])
    ylim([min(lat), max(lat)] + [-5, 5])

    xticklabels(cellfun(@formatlonticks, num2cell(xticks), "UniformOutput", false))
    yticklabels(cellfun(@formatlatticks, num2cell(yticks), "UniformOutput", false))
end
