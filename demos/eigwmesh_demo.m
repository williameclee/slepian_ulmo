function eigwmesh_demo(varargin)
    %% Generating data
    domain = GeoDomain('oceans');
    L = 30;

    [G, V] = glmalpha_new(domain, L, "BeQuiet", true);
    [mesh, lon, lat] = eigwmesh(G, V);
    mesh = mesh / max(mesh(:));

    %% Plotting
    if nargin > 0
        funName = varargin{1};
        figName = ...
            sprintf('eigenvalue-weighted map of Slepian functions'' power (%s)', ...
            upper(funName));
    else
        figName = 'eigenvalue-weighted map of Slepian functions'' power';
    end

    cLim = [0, 1];
    cStep = 0.1;

    figure(999)
    set(gcf, 'Name', figName, 'NumberTitle', 'off')
    clf

    [~, cLevels] = loadcbar(cLim, cStep, 'Normalised power', 'warm');

    hold on
    contourf(lon, lat, mesh, cLevels, 'LineStyle', 'none')
    plotqdm(domain.Lonlat('LonOrigin', 180), 'k')
    hold off
end
