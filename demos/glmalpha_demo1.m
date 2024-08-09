%% GLMALPHA_DEMO1
% This is demo 1 from GLMALPHA.
%
% Authored by:
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12

function glmalpha_demo1
    % Note that using ADDMOUT you can get this back to block-diagonal form
    G = glmalpha;
    [~, ~, ~, bl] = addmout(18);
    imagefnan(G(bl, :));
    axis tight
    difer(G' * G - eye(size(G)))
end