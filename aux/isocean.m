function isOcean = isocean(region)

    if iscell(region)
        region = region{1};
    end

    oceanStrings = ...
        {'oceans', ...
         'atlantic', 'satlantic', 'natlantic', ...
         'pacific', 'spacific', 'npacific', 'indian'};
    isOcean = ismember(region, oceanStrings);
end
