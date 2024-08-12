%% DOMAINNAME
% Returns the full name of a domain given its function name.
%
% Syntax
%   name = domainname(domain)
%   name = domainname(domain, fmt)
%
% Input arguments
%   domain - The domain name
%       The domain name should be the name of a function that returns the
%       domain vertices.
%   fmt - The format of the domain name
%       - 'short' - Short name
%       - 'long' - Long name
%       The default value is 'short'.
%
% Output arguments
%   name - The full name of the domain
%
% Example
%   >>  domainname('namerica')
%   'N America'
%   >>  domainname('namerica', 'long')
%   'North America'
%
% Last modified by
%   2024/08/12, williameclee@arizona.edu (@williameclee)

function domainName = domainname(varargin)
    p = inputParser;
    addRequired(p, 'Domain', ...
        @(x) ischar(x) || isstring(x) || iscell(x));
    addOptional(p, 'Format', 'short', ...
        @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    domain = p.Results.Domain;
    format = p.Results.Format;

    if iscell(domain)
        domain = domain{1};
    end

    domain = lower(domain);

    switch domain
        case 'namerica'

            switch format
                case 'short'
                    domainName = 'N America';
                case 'long'
                    domainName = 'North America';
            end

        case 'samerica'

            switch format
                case 'short'
                    domainName = 'S America';
                case 'long'
                    domainName = 'South America';
            end

        case 'continents'
            domainName = 'All continents';

        case 'oceans'
            domainName = 'All oceans';

        case 'atlantic'

            switch format
                case 'short'
                    domainName = 'Atlantic';
                case 'long'
                    domainName = 'Atlantic Ocean';
            end

        case 'natlantic'

            switch format
                case 'short'
                    domainName = 'N Atlantic';
                case 'long'
                    domainName = 'North Atlantic Ocean';
            end

        case 'satlantic'

            switch format
                case 'short'
                    domainName = 'S Atlantic';
                case 'long'
                    domainName = 'South Atlantic Ocean';
            end

        case 'pacific'

            switch format
                case 'short'
                    domainName = 'Pacific';
                case 'long'
                    domainName = 'Pacific Ocean';
            end

        case 'npacific'

            switch format
                case 'short'
                    domainName = 'N Pacific';
                case 'long'
                    domainName = 'North Pacific Ocean';
            end

        case 'spacific'

            switch format
                case 'short'
                    domainName = 'S Pacific';
                case 'long'
                    domainName = 'South Pacific Ocean';
            end

        case 'indian'

            switch format
                case 'short'
                    domainName = 'Indian';
                case 'long'
                    domainName = 'Indian Ocean';
            end

        otherwise
            domainName = capitalise(domain);

    end

end
