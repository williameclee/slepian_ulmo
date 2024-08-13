%% PLM2XYZ_DEMO3
% This is demo 3 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Last modified by
%   2024/07/25, williameclee@arizona.edu (@williameclee)

function plm2xyz_demo3
% Load some constants - they happen to be the same for both models
a96 = fralmanac('a_EGM96', 'Earth');
GM96 = fralmanac('GM_EGM96', 'Earth');
a2008 = fralmanac('a_EGM2008', 'Earth');
GM2008 = fralmanac('GM_EGM2008', 'Earth');

% Bypass FRALMANAC for efficiency
load(fullfile(getenv('IFILES'), 'EARTHMODELS', 'CONSTANTS', 'SHM')) %#ok<LOAD>

% Load the geopotential coefficients for EGM96
v96 = SHM.EGM96;
% Load the geopotential coefficients for EGM2008
v2008 = SHM.EGM2008_ZeroTide;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP ROW

% Take out the useless error terms and go up to 360 only
v96 = v96(1:addmup(360) - addmup(1), 1:4);
v2008 = v2008(1:addmup(360) - addmup(1), 1:4);

% Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
v96(1, 3) = 0;
v96(:, 3:4) = v96(:, 3:4) * GM96 / a96 .* repmat([(v96(:, 1) - 1) / a96], 1, 2);
v2008(1, 3) = 0;
v2008(:, 3:4) = v2008(:, 3:4) * GM2008 / a2008 .* repmat([(v2008(:, 1) - 1) / a2008], 1, 2);

% Expand, convert to mGal
v96xyz = plm2xyz(v96) * 1e5;
disp(' ')
v2008xyz = plm2xyz(v2008) * 1e5;
disp(' ')

figure(gcf);
clf
[ah, ~] = krijetem(subnum(2, 2));
axes(ah(1));
imagef([], [], v96xyz);
clim([-100 100]);
axis image
cb(1) = colorbar('hor');
% axes(cb(1));
xlabel(cb(1), 'EGM96 free-air anomaly (mgal)')
axes(ah(2));
imagef([], [], v2008xyz);
clim([-100 100]);
axis image
cb(2) = colorbar('hor');
% axes(cb(2));
xlabel(cb(2), 'EGM2008 free-air anomaly (mgal)')
shrink(cb(1:2), 2, 1.5);
movev(cb(1:2), - .075)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM ROW
% Load the geopotential coefficients for EGM96 - again
v96 = SHM.EGM96;
% Load the geopotential coefficients for EGM2008
v2008 = SHM.EGM2008_ZeroTide;
% Take out the useless error terms but increase resolution
v96 = v96(:, 1:4);
v2008 = v2008(1:addmup(1024) - addmup(1), 1:4);

% Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
v96(1, 3) = 0;
v96(:, 3:4) = v96(:, 3:4) * GM96 / a96 .* repmat((v96(:, 1) - 1) / a96, 1, 2);
v2008(1, 3) = 0;
v2008(:, 3:4) = v2008(:, 3:4) * GM2008 / a2008 .* repmat((v2008(:, 1) - 1) / a2008, 1, 2);

% Now specify a certain longitude-latitude grid for partial expansion
c11cmn = [76 14 160 -32];

% Expand, convert to mGal
v96xyz = plm2xyz(v96, [], c11cmn) * 1e5;
disp(' ')
v2008xyz = plm2xyz(v2008, [], c11cmn) * 1e5;
disp(' ')

axes(ah(3));
imagef([], [], v96xyz);
clim([-100 100]);
axis image
cb(3) = colorbar('hor');
% axes(cb(3));
xlabel(cb(3), 'EGM96 free-air anomaly (mgal)')
axes(ah(4));
imagef([], [], v2008xyz);
clim([-100 100]);
axis image
cb(4) = colorbar('hor');
% axes(cb(4));
xlabel(cb(4), 'EGM2008 free-air anomaly (mgal)')
shrink(cb(3:4), 2, 1.5);
movev(cb(3:4), - .075)

fig2print(gcf, 'landscape')
end