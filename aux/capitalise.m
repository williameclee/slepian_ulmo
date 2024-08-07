function capitalisedString = capitalise(inputString)
    if length(inputString) > 1
        capitalisedString = [upper(extractBefore(inputString, 2)), extractAfter(inputString, 1)];
    elseif ~isempty(inputString)
        capitalisedString = upper(inputString);
    else
        capitalisedString = '';
    end

end
