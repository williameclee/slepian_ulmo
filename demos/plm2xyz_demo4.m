%% PLM2XYZ_DEMO4
% This is demo 4 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Last modified by
%   2024/07/25, williameclee@arizona.edu (@williameclee)

function plm2xyz_demo4
    % Bypass FRALMANAC for efficiency
    load(fullfile(getenv('IFILES'), 'EARTHMODELS', 'CONSTANTS', 'SHM')) %#ok<LOAD>

    % Load the topography coefficients for GTM3AR
    vgtm = SHM.GTM3AR;
    % Load the topography coefficients for EGM2008
    v2008 = SHM.EGM2008_Topography;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP ROW

    % Take out the useless error terms and go up to 720 only
    vgtm = vgtm(1:addmup(720), 1:4);
    v2008 = v2008(1:addmup(720), 1:4);

    % Take out the C00-C20 terms to ignore mean and first-order flattening
    vgtm(1:addmup(2), 3:4) = 0;
    v2008(1:addmup(2), 3:4) = 0;

    % Expand, convert to km
    vgtmxyz = plm2xyz(vgtm) * 1e-3;
    v2008xyz = plm2xyz(v2008) * 1e-3;

    figure(gcf);
    clf
    [ah, ~] = krijetem(subnum(2, 2));
    axes(ah(1));
    imagef([], [], vgtmxyz);
    clim([-5 5]);
    axis image
    cb(1) = colorbar('hor');
    % axes(cb(1));
    xlabel(cb(1), 'GTM3AR topography (km)')
    axes(ah(2));
    imagef([], [], v2008xyz);
    clim([-5 5]);
    axis image
    cb(2) = colorbar('hor');
    % axes(cb(2)); 
	xlabel(cb(2), 'EGM2008 topography (km)')
    shrink(cb(1:2), 2, 1.5);
    movev(cb(1:2), - .075)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM ROW
    % Load the topography coefficients for GTM3AR - again
    vgtm = SHM.GTM3AR;
    % Load the topography coefficients for EGM2008
    v2008 = SHM.EGM2008_Topography;
    % Take out the useless error terms but increase resolution
    vgtm = vgtm(:, 1:4);
    v2008 = v2008(1:addmup(1024), 1:4);

    % Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
    vgtm(1:addmup(2), 3:4) = 0;
    v2008(1:addmup(3:4), 3:4) = 0;

    % Now specify a certain longitude-latitude grid for partial expansion
    c11cmn = [76 14 160 -32];

    % Expand, convert to km
    vgtmxyz = plm2xyz(vgtm, [], c11cmn) * 1e-3;
    v2008xyz = plm2xyz(v2008, [], c11cmn) * 1e-3;

    axes(ah(3)); imagef([], [], vgtmxyz); 
	clim([-5 5]); 
	axis image
    cb(3) = colorbar('hor'); 
	% axes(cb(3)); 
	xlabel(cb(3), 'GTM3AR topography (km)')
    axes(ah(4)); 
	imagef([], [], v2008xyz); 
	clim([-5 5]); 
	axis image
    cb(4) = colorbar('hor'); 
	% axes(cb(4)); 
	xlabel(cb(4), 'EGM2008 topography (km)')
    shrink(cb(3:4), 2, 1.5); 
	movev(cb(3:4), - .075)

    fig2print(gcf, 'landscape')
    % I suppose this shows that the two C20 coefficients are very
    % different! Look at Australia.
end
