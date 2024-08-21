function isLeap = isleap(theYear)

    if isdatetime(theYear)
        theYear = theYear.Year;
    end

    theYear = uint16(floor(theYear));
    isLeap = (mod(theYear, 4) == 0) & ((mod(theYear, 100) ~= 0) | (mod(theYear, 400) == 0));

end
