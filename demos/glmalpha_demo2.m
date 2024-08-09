%% GLMALPHA_DEMO2
% This is demo 2 from GLMALPHA.
%
% Authored by:
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12

function glmalpha_demo2
    % Make a demo for Saarimaki et al.
    L = 36;
    [G0, V0, EL0, EM0] = glmalpha(30, L, 1, 0); %#ok<ASGLU>
    [G1, V1, EL1, EM1] = glmalpha(30, L, 1, 1); %#ok<ASGLU>
    [i1, j1] = sort(V1, 'descend'); %#ok<ASGLU>
    [i0, j0] = sort(V0, 'descend'); %#ok<ASGLU>
    N0 = round(sum(V0));
    N1 = round(sum(V1));
    imagefnan(G0(:, 1:N0) * G0(:, 1:N0)')
    imagefnan(G1(:, 1:N1) * G1(:, 1:N1)')

    % Not block sorted
    [G0, V0, EL0, EM0] = glmalpha('africa', L, [], 0); %#ok<ASGLU>
    % Block sorted
    [G1, V1, EL1, EM1] = glmalpha('africa', L, [], 1); %#ok<ASGLU>

    % Fake a coupling kernel a la BCOUPLING
    M = G0(:, 1:N0) * G0(:, 1:N0)';
    M = M .^ 2;

    for l = 0:L
        b = l ^ 2 + 1;
        e = (l + 1) ^ 2;
        % Construct the symmetric part
        for lp = l:L
            bp = lp ^ 2 + 1;
            ep = (lp + 1) ^ 2;
            % The symmetric expression
            ML(l + 1, lp + 1) = sum(sum(M(b:e, bp:ep)));
            ML(lp + 1, l + 1) = ML(l + 1, lp + 1);
        end

    end

    % Forget about the asymmetric scaling and look at a certain cut
    for i = 1:L + 1
        plot(ML(i, :) / ML(i, i))
        j = indeks(find(diff(ML(i, :) / ML(i, i) > 0.05)), 1:2);
        set(gca, 'xtick', sort([i j]), 'xgrid', 'on', ...
            'ytick', [0 0.05 1], 'ygrid', 'on')
        pause
    end

end