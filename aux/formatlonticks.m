function lonL = formatlonticks(lon)
    lon = mod(lon, 360);

    if lon < 180
        lonL = sprintf('%d°E', lon);
    elseif lon > 180
        lonL = sprintf('%d°W', 360 - lon);
    else
        lonL = '180°';
    end

end