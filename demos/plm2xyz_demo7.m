%% PLM2XYZ_DEMO7
% This is demo 7 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Note:
%   This function is not preperly tested. Some inputs and outputs may be missing.
%
% Last modified by
%   2024/07/25, williameclee@arizona.edu (@williameclee)

function plm2xyz_demo7
    % Load the model
    load(fullfile(getenv('IFILES'), ...
        'EARTHMODELS', 'POMME-4', 'pomme-4.2s-nosecular.mat'), 'lmcosi');
    % Bandwidth
    L = 72;

    % Restrict the model to degree L
    lmcosi = lmcosi(1:addmup(L) - addmup(lmcosi(1) - 1), :);

    % Convert to radial-component magnetic field on the reference surface
    lmcosip = plm2mag(lmcosi);

    % Generate a random set of locations inside a spherical patch
    Nd = 1000;
    TH = 15;
    phi0 = 15;
    theta0 = 70;
    [lon, lat] = randpatch(Nd, TH, phi0, theta0);

    % Perform the expansion to the complete grid
    [r1, long, latg] = plm2xyz(lmcosip, 5);

    [LONG, LATG] = meshgrid(long, latg);

    % Now need to convert lmcosi into a format ready for the degree and
    % order ordering implicit in ylm - using ADDMON
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ronm] = addmon(L);

    % Watch for the degree 0 term which isn't in there
    % lmcosi = [0 0 0 0; lmcosi];
    lmcosip = [0 0 0 0; lmcosip];

    % And compare with a full-on gridded version using YLM
    [Y, ~, ~, dems, ~] = ylm([0 L], [], (90 - latg) * pi / 180, long * pi / 180);
    % Don't forget to fix the phase and the normalization
    Y = (Y .* repmat((-1) .^ dems, 1, size(Y, 2))) * sqrt(4 * pi);
    % Expand and reshape
    r2 = reshape(Y' * ...
        lmcosip(2 * size(lmcosip, 1) + ronm), length(latg), length(long));

    % Check again using the irregular version of PLM2XYZ
    r3 = reshape(plm2xyz(lmcosip, LATG(:), LONG(:)), length(latg), length(long));

    % These guys are not identical, but many of them are... and the rest
    % are very close
    difer(r1 - r2, 6)
    difer(r1 - r3, 6)
    difer(r2 - r3, 6)

    tic
    % Now check the expansion on the irregular set
    r4 = plm2xyz(lmcosip, lat, lon);
    toc

    tic
    % And compare with a full-on irregular version using YLM
    [Y, ~, ~, dems, ~] = ...
        ylm([0 L], [], (90 - lat) * pi / 180, lon * pi / 180, [], [], [], 1);
    % Don't forget to fix the phase and the normalization
    Y = (Y .* repmat((-1) .^ dems, 1, size(Y, 2))) * sqrt(4 * pi);

    % Then perform the expansion watching the extra phase factor
    r5 = Y' * lmcosip(2 * size(lmcosip, 1) + ronm);
    toc

    % Check the result
    difer(r4 - r5, 6)

    % Check it out using SCATTER.

    % Plot CRUSTAL model using PLOTPLM, compare with SAARIMAKI1
    lmcosip(1:addmup(15), 3:4) = 0;
    [r1, ~, ~] = plm2xyz(lmcosip, 0.5);
    % Obviously need to check the units here
    imagefnan([0 90], [360 -90], r1, [], [-1000 1000]); axis tight
    plotcont
    minmax(r1 / 1000)
end
