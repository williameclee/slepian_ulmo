%% ENSUREPLMDEGREE
% Pad or truncate a plm matrix to ensure it has the correct degree.
% The input format is assumed to be strictly timefirst (t, lmcosi).
%
% See also
%	SOLVESLE, SOLVEDEGREE1
%
% Notes
%	This is a helper function with limited documentation.%
%
% Authored by
%	2025/05/28, williameclee@arizona.edu (@williameclee)

function plm = ensureplmdegree(plm, L)

    if size(plm, 2) == addmup(L)
        return
    end

    if size(plm, 2) > addmup(L)
        plm = plm(:, 1:addmup(L), :);
        return
    end

    [order, degree] = addmon(L);
    plm(:, 1:addmup(L), 1) = ...
        repmat(reshape(degree, [1, length(degree)]), [size(plm, 1), 1]);
    plm(:, 1:addmup(L), 2) = ...
        repmat(reshape(order, [1, length(degree)]), [size(plm, 1), 1]);
end
