%% CAPITALISE
% Capitalises the first letter of a string.
%
% Syntax
%   sc = capitalise(s)
%
% Input arguments
%   s - The string to capitalise
%
% Output arguments
%   sc - The capitalised string
%
% Last modified by
%   2024/08/12, williameclee@arizona.edu (@williameclee)

function capitalisedString = capitalise(inputString)
    % Empty string
    if isempty(inputString)
        capitalisedString = '';
        return
    end

    % Single character
    if length(inputString) == 1
        capitalisedString = upper(inputString);
        return
    end

    % Multiple characters
    capitalisedString = ...
        [upper(extractBefore(inputString, 2)), ...
         extractAfter(inputString, 1)];

end
