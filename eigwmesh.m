%% EIGWMESH
% Creates a mesh of the eigenvalue-weighted sum of the Slepian functions'
% power on a sphere.
%
% Syntax
%   mesh = eigwmesh(G, V)
%       Finds the mesh of the eigenvalue-weighted sum of the Slepian
%       functions' power.
%   mesh = eigwmesh(G, V, h)
%       Finds the mesh with the specified mesh size.
%   [mesh, lon, lat] = eigwmesh(__)
%       Also returns the longitudes and latitudes of the mesh.
%
% Input arguments
%   G - The localisation matrix of the Slepian functions
%       The matrix is the projection matrix of the Slepian functions onto
%       the spherical harmonics.
%   V - The eigenvalues of the Slepian functions
%   h - The mesh size in degrees
%       The default value is 1.
%
% Output arguments
%   mesh - The mesh of the eigenvalue-weighted sum of the Slepian
%       functions' power
%   lon, lat - The longitudes and latitudes of the mesh
%       The latitudes are ordered from north to south.
%
% Notes
%   This function uses the parallel computing toolbox to speed up the
%   computation.
%
% See also
%   COEF2LMCOSI, PLM2XYZ
%
% Last modified by
%   2024/08/10, williameclee@arizona.edu (@williameclee)

function [mesh, lon, lat] = eigwmesh(varargin)
    %% Initialisation
    % Add path to the auxiliary functions
    addpath(fullfile(fileparts(mfilename('fullpath')), 'aux'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'demos'));

    % Demo
    if isempty(varargin) || strcmpi(varargin{1}, 'demo')
        fprintf('%s running demo for the eigenvalue-weighted power map\n', ...
            upper(mfilename))
        fprintf('Displaying the normalised eigenvalue-weighted power of Slepian functions over the oceans\n')
        equalearth_demo(mfilename)
        return
    end

    % Parse inputs
    [G, V, meshSize] = parseinputs(varargin{:});

    %% Computing the mesh
    lon = 0:meshSize:360;
    lat = 90:-meshSize:-90;
    mesh = zeros([length(lat), length(lon)]);
    % Sum each Slepian function
    parfor i = 1:length(V)
        sphCoeffs = coef2lmcosi(G(:, i), 1);
        meshI = plm2xyz(sphCoeffs, meshSize, "BeQuiet", true);
        mesh = mesh + meshI .^ 2 * V(i);
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addRequired(p, 'G', @(x) ismatrix(x) && isnumeric(x));
    addRequired(p, 'V', @(x) isvector(x) && isnumeric(x));
    addOptional(p, 'MeshSize', 1, @(x) isscalar(x) && isnumeric(x));
    parse(p, varargin{:});
    G = p.Results.G;
    V = p.Results.V;
    meshSize = p.Results.MeshSize;
    % Input validation
    if size(G, 1) ~= size(V, 1)
        error('G and V must have the same number of rows');
    end

    varargout = {G, V, meshSize};
end
