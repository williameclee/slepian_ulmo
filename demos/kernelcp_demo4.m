%% KERNELCP_DEMO4
% This is demo 4 from KERNELCP.
%
% See also
%   KERNELCP, KERNELCP_DEMO1, KERNELCP_DEMO2, KERNELCP_DEMO3, KERNELCP_DEMO5
% Authored by
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12

function kernelcp_demo4
    L = 18;
    [Klmlmp, ~, ~, ~] = kernelc(L, 'antarctica', [], [], 1);
    % Diagonalize
    [C, ~] = eig(Klmlmp);
    % What's up with the eigenvalues? They have small imaginary parts when
    % their magnitude is small.
    difer(eye(size(C)) - (C' * C))
    % If it fails the test, it will be due to this numerical degeneracy!

    % However, the large-eigenvalue ones should still look good!
    % These rotated kernels (whose eigenfunctions should be right on
    % Antarctica without any further ado!) should be orthonormal!
    [dems, dels, ~, ~, mzin] = addmon(L);
    % Preallocate/initialize cell array
    CC = cellnan(J, length(dels), 2);
    % See in LOCALIZATION and GLMALPHA
    for index = 1:size(C)
        % Note that you can undo this using the rinm variable in ADDMON
        CC{index} = reshape(insert(C(:, index), 0, mzin), 2, length(dems))';
        plotplm([dels dems CC{index}], [], [], 2, 0.5); view(145, -65)
        pause
    end

end