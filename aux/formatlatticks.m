function latL = formatlatticks(lat)

    if lat > 0
        latL = sprintf('%d°N', lat);
    elseif lat < 0
        latL = sprintf('%d°S', -lat);
    else
        latL = '0°';
    end

end