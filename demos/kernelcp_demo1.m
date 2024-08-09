%% KERNELCP_DEMO1
% This is demo 1 from KERNELCP.
%
% See also
%   KERNELCP, KERNELCP_DEMO2, KERNELCP_DEMO3, KERNELCP_DEMO4, KERNELCP_DEMO5
% Authored by
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12

function kernelcp_demo1
    L = 12;
    % Construct the kernel for Antarctica rotated to the equator and not
    % rotated back
    K0 = kernelc(L, 'antarctica', [], [], 0);
    % Excessive verification
    xver = 0;
    % Construct the kernel for Antarctica in physical space which means
    % the kernel has been rotated back
    if xver == 1
        % With extra steps
        [K2, ~, K1, K] = kernelc(L, 'antarctica', [], [], 1);
        % This intermediate result better be very nearly the same:
        difer(K - K0)
    else
        % Without the extra steps
        [K2, ~, K1] = kernelc(L, 'antarctica', [], [], 1);
    end

    % Compare the kernels visually
    clf
    strx = 'rank | degree and order';
    stry1 = 'degree | degree and order';
    stry2 = 'degree | degree and order';
    stry3 = 'order | degree and order';
    [ah, ~, ~] = krijetem(subnum(1, 3));

    c11 = [1 1];
    cmn = [(L + 1) ^ 2 (L + 1) ^ 2];

    axes(ah(1))
    % Regular view, degrees indicated
    imagefnan(c11, cmn, K0)
    axis ij
    xlabel(strx)
    % xl(1) = xlabel(strx);
    ylabel(stry1)
    % yl(1) = ylabel(stry1);
    title('The original kernel')
    % t(1) = title('The original kernel');

    axes(ah(2))
    % Regular view, degrees indicated
    % Should look "zonal", almost all on m=0
    imagefnan(c11, cmn, K1); axis ij
    xlabel(strx)
    % xl(2) = xlabel(strx);
    ylabel(stry2)
    % yl(2) = ylabel(stry2);
    title('Once rotated')
    % t(2) = title('Once rotated');

    axes(ah(3))
    % Sort orders and degrees from KERNELC to ADDMOUT
    [~, ~, ~, ~, ~, ~, Km, ~, rinm] = addmon(L);
    % Now [Kl(rinm) Km(rinm)] is like ADDMOUT
    [~, ~, ~, blkm, ~] = addmout(L);
    % Now [Kl(rinm(blkm)),Km(rinm(blkm))] is block sorted
    imagefnan(c11, cmn, K2(rinm(blkm), rinm(blkm))); axis ij
    % blox = Km(rinm(blkm));
    ordshew = [1; find(diff(Km(rinm(blkm)))) + 1];
    % xl(3) = xlabel(strx);
    xlabel(strx)
    % yl(3) = ylabel(stry3);
    ylabel(stry3)
    % t(3) = title('Fully rotated');
    title('Fully rotated')

    % Cosmetics
    longticks(ah)
    set(ah(1:2), 'ytick', addmoff((0:L) - 1) + 1, 'ytickl', 0:L)
    set(ah(3), 'ytick', ordshew, 'ytickl', num2cell(Km(ordshew)))
    set(ah, 'xtick', [1 (L + 1) ^ 2], 'xtickl', [1 (L + 1) ^ 2])
    fig2print(gcf, 'landscape')
    figdisp([], sprintf('antarctica_%i', L))
end
