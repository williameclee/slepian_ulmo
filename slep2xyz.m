%% SLEP2XYZ
% Converts Slepian coefficients to a 2D mesh on a sphere.
%
% Syntax
%   [mesh, lon, lat] = slep2xyz(falpha, domain)
%   [mesh, lon, lat] = slep2xyz(falpha, domain, L)
%   [mesh, lon, lat] = slep2xyz(falpha, domain, L, meshSize)
%   [mesh, lon, lat] = slep2xyz(falpha, domain, L, meshSize, truncation)
%   [mesh, lon, lat, V, N] = slep2xyz(__)
%   slep2xyz('demo')
%
% Input arguments
%   falpha - Slepian expansion coefficients
%   domain - Geographic domain or a latitude-longitude pair
%       - A geographic domain (GeoDomain object).
%       - A string of the domain name.
%       - A cell array of the form {'domain name', buf}.
%       - A N-by-2 matrix of longitude-latitude vertices.
%   L - Bandwidth of the window
%       The default value is the degree of the data.
%   meshSize - The size of the mesh
%       The default value is 1.
%   truncation - The number of Slepian functions to use in the expansion
%       The default value is the Shannon number N.
%
% Output arguments
%   mesh - The 2D mesh on the sphere
%   lon, lat - The longitude and latitude of the mesh
%       The latitude is ordered from north to south
%   V - Eigenvalues of the Slepian functions
%   N - Shannon number
%
% Examples
%   The following example demonstrate the use of SLEP2XYZ (SLEP2PLM + PLM2XYZ):
%   >>  slep2xyz('demo')
%
% Notes
%   The function is essentially a wrapper for SLEP2PLM and PLM2XYZ.
%
% See also
%   PLM2SLEP, SLEP2PLM, PLM2XYZ, XYZ2PLM
%
% Last modified by
%   2024/08/13, williameclee@arizona.edu (@williameclee)

function varargout = slep2xyz(slep, domain, varargin)
    %% Initialisation
    % Add path to the auxiliary functions
    addpath(fullfile(fileparts(mfilename('fullpath')), 'demos'));

    % Demos
    if ischar(slep) || isstring(slep)
        demoId = slep;

        if ~strcmpi(demoId, 'demo')
            error('Unknown demo name ''%s''', demoId);
        end

        slep2plm_demo(mfilename);
        return
    end

    % Parsing inputs
    [slep, domain, L, meshSize, truncation, beQuiet] = ...
        parseinputs(slep, domain, varargin{:});

    %% Computing the mesh
    [lmcosi, V, N] = slep2plm_new(slep, domain, L, "Truncation", truncation);
    [mesh, lon, lat] = plm2xyz(lmcosi, meshSize, "BeQuiet", beQuiet);

    %% Collecting output
    varargout = {mesh, lon, lat, V, N};
end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addRequired(p, 'SlepianCoefficient', @isnumeric);
    addRequired(p, 'Domain', ...
        @(x) ischar(x) || iscell(x) || isa(x, "GeoDomain") || ...
        isnumeric(x));
    addOptional(p, 'L', [], @isnumeric);
    addOptional(p, 'MeshSize', 1, @isnumeric);
    addOptional(p, 'Truncation', [], @isnumeric);
    addParameter(p, 'BeQuiet', false, @islogical);
    parse(p, varargin{:});
    slep = p.Results.SlepianCoefficient;
    domain = p.Results.Domain;
    L = p.Results.L;
    meshSize = p.Results.MeshSize;
    truncation = p.Results.Truncation;
    beQuiet = p.Results.BeQuiet;

    % Convert the domain to a GeoDomain object if applicable
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    % Check if L is valid
    Ldata = sqrt(length(slep)) - 1;

    if L > Ldata
        warning('L is larger than the data supports. Setting L to %d.', Ldata);
    end

    L = min(L, Ldata);

    varargout = {slep, domain, L, meshSize, truncation, beQuiet};
end
