function domainName = domainname(varargin)
    p = inputParser;
    addRequired(p, 'Domain', @(x) ischar(x) || isstring(x) || iscell(x));
    addOptional(p, 'Format', 'short', @(x) ischar(x) || isstring(x));
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
