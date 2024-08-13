%% PLM2XYZ_DEMO2
% This is demo 2 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Last modified by
%   2024/07/25, williameclee@arizona.edu (@williameclee)

function plm2xyz_demo2(varargin)

    if nargin == 0
        meshSize = 0.4;
    else
        meshSize = varargin{1};
    end

    lmax = 10;
    L = 10;
    [m, l, mzero] = addmon(lmax);
    c = randn(addmup(lmax), 2) .* ([l l] .^ (-1));
    c(1) = 3;
    c(mzero, 2) = 0;
    lmcosi = [l m c];
    [r, lon, lat] = plm2xyz(lmcosi, 180 / sqrt(L * (L + 1)));
    tol = length(lon) * length(lat);
    fra = meshSize * 1.5;
    unform = 2;
    [LON, LAT] = meshgrid(lon, lat);

    if unform == 1
        % Not really uniform
        indo = indeks(shuffle((1:tol)'), 1:ceil(fra * tol));
        lonr = LON(indo);
        latr = LAT(indo);
    else
        % Really uniform on the sphere
        [lonr, latr] = randsphere(ceil(fra * tol));
        indo = sub2ind(size(r), ceil(scale(latr, [1 length(lat)])), ...
            ceil(scale(lonr, [1 length(lon)])));
        lonr = LON(indo);
        latr = LAT(indo);
    end

    C6 = xyz2plm(r(indo), L, 'irr', latr, lonr);
    clf
    ah(1) = subplot(221);
    plotplm(r, deg2rad(lon), deg2rad(lat), 1)
    title('Input');
    hold on
    [x, y] = mollweide(deg2rad(lonr), deg2rad(latr));
    plot(x, y, 'o');
    ck = clim;
    cb(1) = colorbar('hor');
    ah(2) = subplot(222);
    rec = plm2xyz(C6);
    plotplm(rec, deg2rad(lon), deg2rad(lat), 1)
    hold on
    plot(x, y, 'o');
    title('Reconstruction');
    clim(ck)
    cb(2) = colorbar('hor');
    ah(3) = subplot(223);
    ps(1) = plot(0:L, plm2spec(lmcosi), 'bo-');
    hold on
    ps(2) = plot(0:L, plm2spec(C6), 'rv-');
    set(ps(1), 'MarkerE', 'b', 'MarkerF', 'r', 'LineW', 2, 'MarkerS', 4)
    set(ps(2), 'MarkerE', 'r', 'MarkerF', 'r', 'LineW', 2, 'MarkerS', 4)
    set(gca, 'Yscale', 'log'); grid
    xlabel('Degree');
    ylabel('Power');
    title(sprintf('Spectral comparison, %s=%3.1f', '\alpha', fra));
    legend('Input', 'Recovered');
    ah(4) = subplot(224);
    difo = abs(r - rec);
    difo(difo < 1e-12) = NaN;
    plotplm(difo, lon * pi / 180, lat * pi / 180, 1)
    hold on
    plot(x, y, 'o');
    % set([p3 p1 p4], 'MarkerE', 'k', 'MarkerF', 'k', 'MarkerS', 2)
    title('Difference > 10^{-12}');
    clim(ck)
    cb(3) = colorbar('hor');
    % set([xl yl], 'FontS', 15)
    longticks(ah, 2)
    movev(cb, - .03)
    kelicol
    fig2print(gcf, 'landscape')
    figdisp('plm2xyz_demo')
end
