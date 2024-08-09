%% KERNELCP_DEMO2
% This is demo 2 from KERNELCP.
%
% See also
%   KERNELCP, KERNELCP_DEMO1, KERNELCP_DEMO3, KERNELCP_DEMO4, KERNELCP_DEMO5
% Authored by
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12

function kernelcp_demo2
    L = 18;
    % Construct the kernel for Antarctica in rotated space
    K0 = kernelc(L, 'antarctica', [], [], 0);

    % Construct the kernel for Antarctica in physical space
    [K2, ~, ~] = kernelc(L, 'antarctica', [], [], 1);

    % Now check the diagonalization of the first kernel
    [C0, V0] = eig(K0);
    [V0, isrt] = sort(sum(real(V0), 1), 'descend');
    C0 = C0(:, isrt(1:addmoff(L)));

    % Now check the diagonalization of the second kernel
    [C2, V2] = eig(K2);
    [V2, isrt] = sort(sum(real(V2), 1), 'descend');
    C2 = C2(:, isrt(1:addmoff(L)));

    % Compare the eigenvalues
    difer(V0 - V2)

    % Compare the eigenfunctions!
    % Clean this up!
    J = 5;
    [dems, dels, ~, ~, mzin] = addmon(L);
    [~, lonc, latc] = antarctica;

    CC = cellnan(J, length(dels), 2);

    parfor index = 1:J
        % This statement reappears exactly in LOCALIZATION
        % And the same arguments appear exactly in KLMLMP2ROT
        CC{index} = plm2rot( ...
            [dels dems reshape(insert(C0(:, index), 0, mzin), 2, length(dems))'], ...
            -lonc, latc, 0);
        % Note that the rotation may lead to some roundoff errors and
        % imaginary parts in the eigenfunctions and eigenvalues!
    end

    % Do the plotting
    clf
    [ah, ha, ~] = krijetem(subnum(J, 3));

    for index = 1:J
        axes(ha(index)) %#ok<LAXES>
        % The non-polar eigenfunctions for Antarctica
        plotslep(C0, index, 2);
        axes(ha(index + J)) %#ok<LAXES>
        % The near-polar eigenfunctions for Antarctica
        p = plotslep(C2, index, 2);
        axes(ha(index + 2 * J)) %#ok<LAXES>
        % The rotated eigenfunctions of the non-polar kernel
        d = plotplm(CC{index}, [], [], 4, 1);
        set(ah, 'clim', halverange(d, 100))
        % Let us appraise the comparison also - knowing that the sign remains
        % arbitrary
        a = minmax(abs(d - p) ./ d);
        b = minmax(abs(d + p) ./ d);
        fprintf('The min/max relative difference is %8.3e | %8.3e', ...
            [max([a(1) b(1)]) min([a(2) b(2)])])
    end

    % Cosmetics
    nolabels(ha(J + 1:end), 2)
    nolabels(ah(1:end - 3), 1)
    longticks(ah)
    fig2print(gcf, 'tall')
end