%% ISOCEAN
% Determine if a domain is an ocean basin.
%
% Syntax
%   flag = isocean(domain)
%
% Input arguments
%   domain - The domain to check
%       The domain can be a string, a cell array of strings, or a GeoDomain
%       object.
%
% Output arguments
%   flag - A flag indicating whether the domain is an ocean basin
%
% Notes
%   The following domains are considered ocean basins:
%       oceans, arctic, indian, pacific, spacific, npacific, atlantic,
%       satlantic, natlantic
%
% Last modified by
%   2024/08/10, williameclee@arizona.edu (@williameclee)

function isOcean = isocean(domain)

    if iscell(domain)
        domain = domain{1};
    elseif isa(domain, 'GeoDomain')
        domain = domain.Domain;
    end

    oceanStrings = ...
        {'oceans', 'arctic', 'indian', ...
         'pacific', 'spacific', 'npacific', ...
         'atlantic', 'satlantic', 'natlantic'};
    isOcean = ismember(domain, oceanStrings);
end
