function varargout = domainspecs(varargin)
    upscaleDefault = [];
    bufDefault = 0;
    inclangDefault = 90;
    moreBufferDefault = {};
    p = inputParser;
    addRequired(p, 'Domain', @ischar);
    addOptional(p, 'Upscale', upscaleDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || isempty(x));
    addOptional(p, 'Buffer', bufDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || isempty(x));
    addOptional(p, 'Inclang', inclangDefault, ...
        @(x) (isnumeric(x) && length(x) <= 2) || isempty(x));
    addOptional(p, 'MoreBuffer', moreBufferDefault, ...
        @(x) (iscell(x) && mod(length(x), 2) == 0) || isempty(x));

    if iscell(varargin{1})
        varargin = [varargin{1}, varargin(2:end)];
    end

    parse(p, varargin{:});

    domain = p.Results.Domain;
    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    inclang = p.Results.Inclang;
    moreBuf = p.Results.MoreBuffer;

    if isempty(upscale)
        upscale = upscaleDefault;
    end

    if isempty(buf)
        buf = bufDefault;
    end

    if isempty(inclang)
        inclang = inclangDefault;
    end

    if isscalar(inclang)
        inclang = [-inclang, inclang];
    end

    if isempty(moreBuf)
        moreBuf = moreBufferDefault;
    end

    moreBufDomain = moreBuf(1:2:end);
    moreBufBuf = moreBuf(2:2:end);

    if any(cellfun(@(x) ~ischar(x) || isstring(x), moreBufDomain))
        error('Domain names must be strings');
    else
        moreBufDomain = cellfun(@(x) char(lower(x)), moreBufDomain, ...
            'UniformOutput', false);
    end

    if any(cellfun(@(x) ~(isnumeric(x) && isscalar(x)), moreBufBuf))
        error('Buffer values must be scalars');
    elseif any(cellfun(@(x) x < 0, moreBufBuf))
        error('Buffer values must be non-negative');
    end

    % Remove zero buffer values
    isZeroBuf = cellfun(@(x) x == 0, moreBufBuf);
    moreBufDomain = moreBufDomain(~isZeroBuf);
    moreBufBuf = moreBufBuf(~isZeroBuf);

    moreBufDomain = moreBufDomain(end:-1:1);
    moreBufBuf = moreBufBuf(end:-1:1);
    
    % Remove duplicates
    [moreBufDomain, idx] = unique(moreBufDomain);
    moreBufBuf = moreBufBuf(idx);

    % Sort by domain name
    [~, idx] = sort(moreBufDomain);
    moreBufDomain = moreBufDomain(idx);
    moreBufBuf = moreBufBuf(idx);

    moreBuf = [moreBufDomain(:)'; moreBufBuf(:)'];
    moreBuf = moreBuf(:)';

    domainSpecs = ...
        {domain, 'Upscale', upscale, 'Buffer', buf, ...
         'Inclang', inclang, 'MoreBuffer', moreBuf};

    varargout = {domainSpecs, upscale, buf, inclang, moreBuf};
end
