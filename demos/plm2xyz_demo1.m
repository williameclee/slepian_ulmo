%% PLM2XYZ_DEMO1
% This is demo 1 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Notes
%   This function is not preperly tested. Some inputs and outputs may be missing.
%
% Last modified by
%   2024/07/25, williameclee@arizona.edu (@williameclee)

function plm2xyz_demo1
    lmax = 30; 
	[m, l, mzero] = addmon(lmax);
    c = randn(addmup(lmax), 2) .* ([l l] .^ (-1));
    c(1) = 3; 
	c(mzero, 2) = 0; 
	lmcosi = [l m c];
    L = 30;
    [r, ~, ~] = plm2xyz(lmcosi, 180 / sqrt(L * (L + 1)));
    C1 = xyz2plm(r, L, 'simpson');
    C2 = xyz2plm(r, L, 'gl');
    C3 = xyz2plm(r, L, 'im');
    lat = linspace(90, -90, size(r, 1));
    C4 = xyz2plm(r, L, 'im', lat);
    C5 = xyz2plm(r, L, 'irr');
    clf
    ah(1) = subplot(211);
    p1(1) = plot(abs(lmcosi(:, 3) - C1(1:addmup(lmax), 3)), 'b+-'); hold on
    p1(2) = plot(abs(lmcosi(:, 3) - C2(1:addmup(lmax), 3)), 'rv-');
    p1(3) = plot(abs(lmcosi(:, 3) - C3(1:addmup(lmax), 3)), 'kx-');
    p1(4) = plot(abs(lmcosi(:, 3) - C4(1:addmup(lmax), 3)), 'go-');
    p1(5) = plot(abs(lmcosi(:, 3) - C5(1:addmup(lmax), 3)), 'mo-'); hold off
    ylim([-0.1 1] * 1e-14)
    yl(1) = ylabel('Absolute error');
    xl(1) = xlabel('Cumulative degree and order');
    legend('simpson', 'gl', 'im', 'imlat', 'irr')

    ah(2) = subplot(212);
    p2(1) = plot(abs(lmcosi(:, 3) - C1(1:addmup(lmax), 3)), 'b+-'); hold on
    p2(2) = plot(abs(lmcosi(:, 3) - C2(1:addmup(lmax), 3)), 'rv-');
    p2(3) = plot(abs(lmcosi(:, 3) - C3(1:addmup(lmax), 3)), 'kx-');
    p2(4) = plot(abs(lmcosi(:, 3) - C4(1:addmup(lmax), 3)), 'go-');
    p2(5) = plot(abs(lmcosi(:, 3) - C5(1:addmup(lmax), 3)), 'mo-'); hold off
    axis tight
    yl(2) = ylabel('Absolute error');
    xl(2) = xlabel('Cumulative degree and order');

    longticks(ah, 2)
    % set([p1 p2], 'MarkerS', 4)
    % set([xl yl], 'FontS', 15)
    fig2print(gcf, 'landscape')
    figdisp
end
