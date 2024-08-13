%% GLMALPHA_DEMO1
% This is demo 1 from GLMALPHA.
%
% Last modified by
%   2024/07/12, williameclee@arizona.edu (@williameclee)

function glmalpha_demo1
    % Note that using ADDMOUT you can get this back to block-diagonal form
    G = glmalpha;
    [~, ~, ~, bl] = addmout(18);
    imagefnan(G(bl, :));
    axis tight
    difer(G' * G - eye(size(G)))
end