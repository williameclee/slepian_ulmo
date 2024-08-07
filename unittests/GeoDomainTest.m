% Test 1: Basic functionality
domain = GeoDomain('oceans');
assert(strcmp(domain.Domain, 'oceans'));
assert((domain.Upscale) == 0);
assert((domain.Buffer) == 0);
assert(isequal(domain.Latlim, [-90, 90]));
assert(isempty(domain.MoreBuffers));
assert(~domain.NearBy);

% Test 2: Assigning default parameters
domain = GeoDomain('indian', 'DefaultParams', true);
assert(strcmp(domain.Domain, 'indian'));
assert((domain.Upscale) == 0);
assert((domain.Buffer) == 0);
assert(isequal(domain.Latlim, [-60, 60]));
assert(isequal(domain.MoreBuffers, {'earthquakes', 10}));
assert(~domain.NearBy);

% Test 3: Assigning custom parameters
domain = GeoDomain('invaliddomain', 'Upscale', -1, 'Buffer', 1, ...
    'Latlim', [-30, 30], 'MoreBuffers', {'earthquakes', 10}, ...
    'NearBy', true);
assert(strcmp(domain.Domain, 'invaliddomain'));
assert((domain.Upscale) == 0);
assert((domain.Buffer) == 1);
assert(isequal(domain.Latlim, [-30, 30]));
assert(isequal(domain.MoreBuffers, {'earthquakes', 10}));
assert(~domain.NearBy);

% Test 4: Multiple domains
domain = GeoDomain(["oceans", "indian"], 'DefaultParams', true);
assert(strcmp(domain(1).Domain, 'oceans'));
assert(strcmp(domain(2).Domain, 'indian'));
assert((domain(1).Upscale) == 0);
assert((domain(2).Upscale) == 0);
assert((domain(1).Buffer) == 0);
assert((domain(2).Buffer) == 0);
assert(isequal(domain(1).Latlim, [-90, 90]));
assert(isequal(domain(2).Latlim, [-60, 60]));
assert(isequal(domain(1).MoreBuffers, {'earthquakes', 10}));
assert(isequal(domain(2).MoreBuffers, {'earthquakes', 10}));
assert(~domain(1).NearBy);
assert(~domain(2).NearBy);

domain = GeoDomain("oceans", 'DefaultParams', true, "Buffer", 1:6);
assert(all(strcmp([domain.Domain], 'oceans')));
assert(all([domain.Upscale] == 0));
assert(isequal([domain.Buffer], 1:6));