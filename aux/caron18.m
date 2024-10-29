dataFolder = fullfile(getenv('IFILES'), 'GIA', 'Caron18');
inputFile = 'GIA_maps_Caron_et_al_2018.txt';
inputPath = fullfile(dataFolder, inputFile);

gia = readmatrix(inputPath, "HeaderLines", 6);
lonR = gia(:, 2);
latR = 90 - gia(:, 1);
giaVlmR = gia(:, 3);
giaVlmErrR = gia(:, 4);
giaGeoidR = gia(:, 5);
giaGeoidErrR = gia(:, 6);

lon = 0:1:360;
lat = -90:1:90;
[lonn, latt] = meshgrid(lon, lat);
F = scatteredInterpolant(lonR, latR, giaVlmR, 'linear', 'linear');
giaVlmLonlat = F(lonn, latt);
F = scatteredInterpolant(lonR, latR, giaVlmErrR, 'linear', 'linear');
giaVlmErrLonlat = F(lonn, latt);
F = scatteredInterpolant(lonR, latR, giaGeoidR, 'linear', 'linear');
giaGeoidLonlat = F(lonn, latt) / 1000; % mm/yr -> m/yr
F = scatteredInterpolant(lonR, latR, giaGeoidErrR, 'linear', 'linear');
giaGeoidErrLonlat = F(lonn, latt) / 1000; % mm/yr -> m/yr

giaVlmPlm = xyz2plm(flip(giaVlmLonlat), 60);
giaVlmErrPlm = xyz2plm(flip(giaVlmErrLonlat), 60);
giaGeoidPlm = xyz2plm(flip(giaGeoidLonlat), 60);
giaSdPlm = plm2pot(giaGeoidPlm, [], [], [], 4);
giaSdLonlat = flip(plm2xyz(giaSdPlm, 1));
giaGeoidErrPlm = xyz2plm(flip(giaGeoidErrLonlat), 60);
giaSdErrPlm = plm2pot(giaGeoidErrPlm, [], [], [], 4);
giaSdErrLonlat = flip(plm2xyz(giaSdErrPlm, 1));

%% Save data
% Surface mass density
lmcosiM = giaSdPlm;
lmcosiU = giaSdPlm;
lmcosiU(:, 3:4) = giaSdPlm(:, 3:4) + giaSdErrPlm(:, 3:4);
lmcosiL = giaSdPlm;
lmcosiL(:, 3:4) = giaSdPlm(:, 3:4) - giaSdErrPlm(:, 3:4);

save(fullfile(dataFolder, 'Caron18_SD.mat'), 'lmcosiM', 'lmcosiU', 'lmcosiL')

% Vertical land motion
lmcosiM = giaVlmPlm;
lmcosiU = giaVlmPlm;
lmcosiU(:, 3:4) = giaVlmPlm(:, 3:4) + giaVlmErrPlm(:, 3:4);
lmcosiL = giaVlmPlm;
lmcosiL(:, 3:4) = giaVlmPlm(:, 3:4) - giaVlmErrPlm(:, 3:4);

save(fullfile(dataFolder, 'Caron18_VLM.mat'), 'lmcosiM', 'lmcosiU', 'lmcosiL')

%% Plotting
figure(1)
clf

contourf(lon, lat, giaGeoidLonlat, 100, 'LineStyle', 'none')
colorbar

figure(2)
clf

contourf(lon, lat, giaSdLonlat, 100, 'LineStyle', 'none')
colorbar
