%% KERNELORDER
%  Returns the index of lmcosi corresponding to the kernel of order L.
%
% See also
%   KERNELCP_NEW
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)

function j = kernelorder(L)
    do = ones([(L + 1) * (L + 2) / 2, 2]);
    do(((0:L) .* (1:L + 1)) / 2 + 1, 2) = 0;
    do = do';
    do(do ~= 0) = 1:nnz(do(:) ~= 0);
    do = do';
    do = do(:);
    [~, j] = sort(do);
    j = j(L + 2:end);
end
