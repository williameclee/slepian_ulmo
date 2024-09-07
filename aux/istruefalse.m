%% ISTRUEFALSE
% Last modified by
%   2024/09/06, williameclee@arizona.edu (@williameclee)

function val = istruefalse(arg, varargin)
    p = inputParser;
    addOptional(p, 'Auto', false, @islogical);
    parse(p, varargin{:});

    useAuto = p.Results.Auto;

    if islogical(arg)
        val = true;
        return
    end

    if isnumeric(arg)
        % The value of the argument must be a integer
        if floor(arg) == arg
            val = true;
        else
            val = false;
        end

        return
    end

    if ischar(arg) || isstring(arg)
        arg = lower(char(arg));

        if ismember(arg, {'true', 'false', 'on', 'off'})
            val = true;
            return
        end

        if useAuto && strcmp(arg, 'auto')
            val = true;
            return
        end

    end

    val = false;

end
