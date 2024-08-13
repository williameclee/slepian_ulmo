%% CONDDEFVAL
% Conditionally defines a value based on whether it is empty.
%
% Syntax
%   v = conddefval(v, v0)
%
% Input arguments
%   v - The input value to be checked
%   v0 - The default value to be returned if the input value is empty
%
% Output argument
%   v - The resulting value based on the condition.
%
% Example
%   The following example will return v = 10:
%   >>  v = [];
%   >>  v = conddefval(v, 10);
%
% Last modified by
%   2024/08/13, williameclee@arizona.edu (@williameclee)

function v = conddefval(v, v0, varargin)
    p = inputParser;
    addRequired(p, 'v');
    addRequired(p, 'v0');
    addParameter(p, 'IncludeNan', false, @(x) islogical(x) || isnumeric(x));
    parse(p, v, v0, varargin{:});
    v = p.Results.v;
    v0 = p.Results.v0;
    includeNan = logical(p.Results.IncludeNan);

    if ~isempty(v) && ~(includeNan && isnumeric(x) && isnan(v))
        return
    end

    v = v0;

end
