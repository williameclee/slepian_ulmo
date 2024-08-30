[dataPlmt, date] = mascon2plmt;
[mesh, lon, lat] = plm2xyz(dataPlmt(:, :, 1));

contourf(lon, lat, mesh, 100, 'LineStyle', 'none');