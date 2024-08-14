%% SLEP2PLM_DEMO
% Is a demo for the functions SLEP2PLM and PLM2SLEP
% It localises a surface mass density anomaly field (from the first month
% of GRACE data) to the North Pacific Ocean using the Slepian basis
% functions.
%
% See also
%   SLEP2PLM, PLM2SLEP
%
% Last modified by
%   2024/08/13, williameclee@arizona.edu (@williameclee)

function slep2plm_demo(varargin)
    %% Initialisation
    % Parse inputs
    p = inputParser;
    addOptional(p, 'FunctionName', 'slep2plm', ...
        @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    funName = p.Results.FunctionName;

    % Load the sample data
    domain = GeoDomain('npacific', 'Buffer', 1, 'DefaultParams', true);
    L = 60;
    meshSize = 1;

    load('sampledata.mat', 'Plm');
    Plm(:, 3:4) = Plm(:, 3:4) / 1e4; % kg/m^2 -> 10^4kg/m^2

    %% Calculation
    % Project the data to the Slepian basis
    slep = plm2slep_new(Plm, domain, L);
    % Convert the spherical harmonic coefficients to meshes
    [mesh, lon, lat] = plm2xyz(Plm, meshSize, "BeQuiet", true);

    % Localise the with Slepian functions and convert them to meshes
    meshLcl = slep2xyz(slep, domain, L, meshSize, ...
        "BeQuiet", true);

    %% Plotting
    % Load the coastline
    coastLonlat = gshhscoastline('l');

    % Colurmap settings
    cLim = [-150, 200];
    cStep = 50/3;
    cLevels = cLim(1):cStep:cLim(2);

    % Miscelaneous settings
    figName = sprintf('Demo of %s', upper(funName));
    cbarTitle = 'Surface mass density anomaly [10^4 kg/m^2]';

    figure(999)
    set(gcf, "Name", figName, "NumberTitle", 'off')
    clf

    % Map before localisation
    subplot(1, 2, 1)
    loadbasemap(domain, "PrjMethod", 'eqaazim')
    loadcbar(cLevels, "Title", cbarTitle, "Location", 'southoutside')
    title('Before localisation')

    contourfm(lat, lon, mesh, cLevels, "LineStyle", 'none')
    plotm(coastLonlat(:, 2), coastLonlat(:, 1), 'k', "LineWidth", 0.5)
    plotm(domain.Lat("Anchors", true), domain.Lon("Anchors", true), ...
        'k', "LineWidth", 1.5)

    % Map after localisation
    subplot(1, 2, 2)
    loadbasemap(domain, "PrjMethod", 'eqaazim')
    loadcbar(cLevels, "Title", cbarTitle, "Location", 'southoutside')
    title('After localisation')

    contourfm(lat, lon, meshLcl, cLevels, "LineStyle", 'none')
    plotm(coastLonlat(:, 2), coastLonlat(:, 1), 'k', "LineWidth", 0.5)
    plotm(domain.Lat("Anchors", true), domain.Lon("Anchors", true), ...
        'k', "LineWidth", 1.5)
end
