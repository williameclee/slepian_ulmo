function time = year2date(timeInYears)
    theYear = floor(timeInYears);
    timeInDays = (timeInYears - theYear) .* ...
        (365 + double(isleap(theYear)));
    theDay = floor(timeInDays);
    timeInHours = (timeInDays - theDay) * 24;
    theHour = floor(timeInHours);
    timeInMinutes = (timeInHours - theHour) * 60;
    theMinute = floor(timeInMinutes);
    theSecond = round((timeInMinutes - theMinute) * 60);

    time = datetime(theYear, 1, theDay, theHour, theMinute, theSecond);
end