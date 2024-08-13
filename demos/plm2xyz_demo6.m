%% PLM2XYZ_DEMO6
% This is demo 6 from PLM2XYZ.
%
% See also
%   PLM2XYZ
%
% Note
%   This function is not preperly tested. Some inputs and outputs may be missing.
%
% Last modified by
%   2024/07/25, williameclee@arizona.edu (@williameclee)

function Outputs = plm2xyz_demo6(nOutputs)
    v = fralmanac('EGM2008_Topography', 'SHM');
    v = v(1:addmup(1000) - addmup(v(1) - 1), :);
    [r, lon, lat, Plm, h] = plm2xyz(v, [], [8 45 20 37]);

    if nOutputs == 0
        [~, ~, ~] = plotplm(r, lon * pi / 180, lat * pi / 180, 4);
        disp('Adjust the axes!!')
    end

    % Now with Dongs' which looks good, passed all tests 4/21/2010
    % [rdw,londw,latdw]=plm2xyzdw(v,[],[45 37 8 20],0,0);

    % Prepare output
    Outputs = {r, lon, lat, Plm, h};
end
