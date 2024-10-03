domainList = {'oceans', 'npacific', 'spacific', 'natlantic', 'satlantic', 'indian'};
Altimetry = ludwigsenaltimetry(domainList);
GraceLudwigsen = ludwigsengrace(domainList);
timeRange = ...
    [min([Altimetry.Time(1), GraceLudwigsen.Time(1)]), ...
          max([Altimetry.Time(end), GraceLudwigsen.Time(end)])];
GraceOur = ourgrace(domainList, timeRange);

for iDomain = 1:length(domainList)
    domain2plot = domainList{iDomain};
    GraceOur.(capitalise(domain2plot)) = ...
        GraceOur.(capitalise(domain2plot)) - mean(GraceOur.(capitalise(domain2plot)), "omitmissing") ...
        + mean(GraceLudwigsen.(capitalise(domain2plot)), "omitmissing");

    figure(1)
    clf
    hold on
    plot(Altimetry.Time, Altimetry.([capitalise(domain2plot), '_StrAlt']), ...
        "Color", 'k', ...
        "DisplayName", "Altimetry (Ludwigsen et al., 2024)");
    plot(GraceLudwigsen.Time, GraceLudwigsen.(capitalise(domain2plot)), ...
        "Color", cc('b2'), "Linewidth", 1, ...
        "DisplayName", "GRACE (Ludwigsen et al., 2024)");
    plotwithgap(GraceOur.Time, GraceOur.(capitalise(domain2plot)), [], ...
        "Color", cc('r2'), "Linewidth", 1, ...
        "DisplayName", "GRACE (This study)");
    hold off

    title(domainname(domain2plot, 'long'))
    ylabel('Ocean mass change [mm]')
    legend("Box", 'off', "Location", 'southeast')

    ylim([-20, 60])

    set(gca, "Box", 'on')
    ann = annotation('rectangle', [0, 0, 1, 1], ...
        "FaceColor", 'none', "EdgeColor", [1, 1, 1] * 0.99);
    exportgraphics(gcf, fullfile(getenv('FIGURES'), '2024', ...
        ['grace_comparison_', domain2plot, '.pdf']), "ContentType", 'vector', "BackgroundColor", 'none');
    delete(ann)
end

function Altimetry = ludwigsenaltimetry(domainList)
    dataPathO = mfilename('fullpath');
    dataPathO = [dataPathO, '.mat'];

    if exist(dataPathO, 'file')
        load(dataPathO, 'Altimetry');

        if exist('Altimetry', 'var')
            return
        end

    end

    dataPath = fullfile(getenv('IFILES'), 'OTHERS', 'Ludwigsen2024.xlsx');

    nDomains = length(domainList);
    domainSheetName = {'Global', 'North Pacific', 'South Pacific', 'North Atlantic', 'South Atlantic', 'Indian Ocean'};

    readtable(dataPath, "Sheet", domainSheetName{1}).time(:);
    Altimetry.Time = year2date(readtable(dataPath, "Sheet", domainSheetName{1}).time);

    for iDomain = 1:nDomains
        AltimetryBasin = readtable(dataPath, "Sheet", domainSheetName{iDomain});
        Altimetry.([capitalise(domainList{iDomain}), '_StrAlt']) = AltimetryBasin.steAlt;
        Altimetry.([capitalise(domainList{iDomain}), '_StrAlt_NGia']) = AltimetryBasin.steAlt - AltimetryBasin.GIA;
        Altimetry.([capitalise(domainList{iDomain}), '_GIA']) = AltimetryBasin.GIA;
    end

    try
        save(dataPathO, 'Altimetry', '-append');
    catch
        save(dataPathO, 'Altimetry');
    end

end

function LudwigsenGrace = ludwigsengrace(domainList)
    dataPathO = mfilename('fullpath');
    dataPathO = [dataPathO, '.mat'];

    if exist(dataPathO, 'file')
        load(dataPathO, 'LudwigsenGrace');

        if exist('LudwigsenGrace', 'var')
            return
        end

    end

    dataPath = fullfile(getenv('IFILES'), 'OTHERS', 'Ludwigsen2024.xlsx');

    nDomains = length(domainList);
    domainSheetName = {'Global', 'North Pacific', 'South Pacific', 'North Atlantic', 'South Atlantic', 'Indian Ocean'};

    readtable(dataPath, "Sheet", domainSheetName{1}).time(:);
    LudwigsenGrace.Time = year2date(readtable(dataPath, "Sheet", domainSheetName{1}).time);

    for iDomain = 1:nDomains
        AltimetryBasin = readtable(dataPath, "Sheet", domainSheetName{iDomain});
        LudwigsenGrace.(capitalise(domainList{iDomain})) = AltimetryBasin.GRACE;
    end

    try
        save(dataPathO, 'LudwigsenGrace', '-append');
    catch
        save(dataPathO, 'LudwigsenGrace');
    end

end

function Grace = ourgrace(domainList, timeRange)
    dataPathO = mfilename('fullpath');
    dataPathO = [dataPathO, '.mat'];

    if exist(dataPathO, 'file')
        load(dataPathO, 'Grace');

        if exist('Grace', 'var')
            return
        end

    end

    for i = 1:length(domainList)
        domain = GeoDomain(domainList{i}, "Buffer", 1, "DefaultParams", true);
        [date, total] = grace2trend(domain, 60, timeRange, [3, 365.0, 181.0], {'CSR', 'RL06', 60}, "GIACorrection", 'mascon');

        Grace.Time = date;
        Grace.(capitalise(domainList{i})) = mass2weq(total, domain, 'seawater', 'mm');
    end

    try
        save(dataPathO, 'Grace', '-append');
    catch
        save(dataPathO, 'Grace');
    end

end
