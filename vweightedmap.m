function [mesh, lon, lat] = vweightedmap(G, V, varargin)
    p = inputParser;
    addRequired(p, 'G', @(x) ismatrix(x) && isnumeric(x));
    addRequired(p, 'V', @(x) isvector(x) && isnumeric(x));
    addOptional(p, 'MeshSize', 1, @(x) isscalar(x) && isnumeric(x));
    parse(p, G, V, varargin{:});
    G = p.Results.G;
    V = p.Results.V;
    meshSize = p.Results.MeshSize;

    if size(G, 1) ~= size(V, 1)
        error('G and V must have the same number of rows');
    end

    lon = 0:meshSize:360;
    lat = 90:-meshSize:-90;
    mesh = zeros([length(lat), length(lon)]);

    parfor i = 1:length(V)
        sphCoeffs = coef2lmcosi(G(:, i), 1);
        meshI = plm2xyz(sphCoeffs, meshSize, "BeQuiet", true);
        mesh = mesh + meshI .^ 2 * V(i);
    end

end
