%% PLM2XYZ_DEMO5
% This is demo 5 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Authored by
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-25

function plm2xyz_demo5
    % Checks that the latest modifications have been put in right
    % Load some constants - they happen to be the same for both models
    a96 = fralmanac('a_EGM96', 'Earth');
    GM96 = fralmanac('GM_EGM96', 'Earth');

    % Bypass FRALMANAC for efficiency
    load(fullfile(getenv('IFILES'), 'EARTHMODELS', 'CONSTANTS', 'SHM')) %#ok<LOAD>

    % Load the geopotential coefficients for EGM96
    v96 = SHM.EGM96;

    % Take out the useless error terms but increase resolution
    v96 = v96(:, 1:4);

    % Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
    v96(1, 3) = 0;
    v96(:, 3:4) = v96(:, 3:4) * GM96 / a96 .* repmat((v96(:, 1) - 1) / a96, 1, 2);

    % In chunks of 72
    r1 = plm2xyz(v96, [], [], 90);
    % In chunks of 144
    r2 = plm2xyz(v96, [], [], 180);
    % All at once
    r3 = plm2xyz(v96, [], [], 360);

    % Check that the results are identical
    difer(r1 - r2); difer(r2 - r3); difer(r3 - r1)
end
