%% Getting GIA displacement trends
L = 60;
timeRange = datetime([2003, 2022], [1, 12], [1, 31]);

domainList = {'oceans', 'npacific', 'spacific', 'natlantic', 'satlantic', 'indian'};
nDomains = length(domainList);

[modelList, nModels] = getmodellist;

GiaCorrAltimetry = getdisplacementdata(domainList, modelList, L);
AltimetryRaw = getaltimetrydata(domainList);
GiaCorrAltimetry(end + 1, :) = AltimetryRaw(2, :);
Altimetry = AltimetryRaw{1, :} - GiaCorrAltimetry;

clear AltimetryRaw

modelList{end + 1} = 'ice6gd_vm5a';
nModels = nModels + 1;

%% Getting GRACE trends
[GraceRaw, GiaCorrMass] = getgracedata(domainList, modelList, L, timeRange);
Grace = GraceRaw{1, :} - GiaCorrMass;
clear GraceRaw

%% Plotting
domainDisplaynames = ...
    cellfun(@domainname, domainList, 'UniformOutput', false);

figure(1)
clf

plottotalcomparison(Altimetry, Grace, modelList, domainDisplaynames)

exportgraphics(gcf, fullfile(getenv('FIGURES'), 'altimetry_grace_comparison.png'))

figure(2)
clf

plotgiacomparison(GiaCorrAltimetry, modelList, domainDisplaynames, 'o')

exportgraphics(gcf, fullfile(getenv('FIGURES'), 'altimetry_gia_comparison.png'))

clf
plotgiacomparison(GiaCorrMass, modelList, domainDisplaynames, '+')

exportgraphics(gcf, fullfile(getenv('FIGURES'), 'grace_gia_comparison.png'))

figure(3)
clf

hold on

for iModel = 1:nModels
    model = erase(modelList{iModel}, {'Steffen_', '-ice'});
    model = strsplit(model, '_');

    if contains(model{1}, 'ice6gd')
        mcolor = cc('g3');
    elseif contains(model{1}, 'ice6g')
        mcolor = cc('b3');
    elseif contains(model{1}, 'ice7g')
        mcolor = cc('r3');
    elseif contains(model{1}, 'anu')
        mcolor = cc('y3');
    else
        mcolor = 'y';
    end

    plot(1:nDomains, Grace{iModel, :} - Altimetry{iModel, :}, 'o-', "Color", mcolor, ...
        "HandleVisibility", 'off')
end

plot(nan, nan, 'o-', "Color", cc('g3'), "DisplayName", 'ice6g\_d')
plot(nan, nan, 'o-', "Color", cc('b3'), "DisplayName", 'ice6g\_c')
plot(nan, nan, 'o-', "Color", cc('r3'), "DisplayName", 'ice7g')
plot(nan, nan, 'o-', "Color", cc('y3'), "DisplayName", 'anu')
hold off

yline(0, 'k--', "HandleVisibility", 'off')

ylabel('Ocean mass trend (GRACE - Altimetry) [mm/yr]')
legend("Box", 'off')

xticks(1:nDomains)
xticklabels(domainDisplaynames)
xlim([1, nDomains] + [-0.5, 0.5])
box on

exportgraphics(gcf, fullfile(getenv('FIGURES'), 'altimetry_grace_comparison_diff.png'))

figure(4)
clf

hold on

for iModel = 1:nModels
    model = erase(modelList{iModel}, {'Steffen_', '-ice'});
    model = strsplit(model, '_');

    if contains(model{1}, 'ice6gd')
        mcolor = cc('g3');
    elseif contains(model{1}, 'ice6g')
        mcolor = cc('b3');
    elseif contains(model{1}, 'ice7g')
        mcolor = cc('r3');
    elseif contains(model{1}, 'anu')
        mcolor = cc('y3');
    else
        mcolor = 'y';
    end

    plot(1:nDomains, abs(Grace{iModel, :} - Altimetry{iModel, :}), 'o-', "Color", mcolor, ...
        "HandleVisibility", 'off')
end

plot(nan, nan, 'o-', "Color", cc('g3'), "DisplayName", 'ice6g\_d')
plot(nan, nan, 'o-', "Color", cc('b3'), "DisplayName", 'ice6g\_c')
plot(nan, nan, 'o-', "Color", cc('r3'), "DisplayName", 'ice7g')
plot(nan, nan, 'o-', "Color", cc('y3'), "DisplayName", 'anu')
hold off

ylabel('Ocean mass trend (abs GRACE - Altimetry) [mm/yr]')
legend("Box", 'off')

xticks(1:nDomains)
xticklabels(domainDisplaynames)
xlim([1, nDomains] + [-0.5, 0.5])
box on

exportgraphics(gcf, fullfile(getenv('FIGURES'), 'altimetry_grace_comparison_absdiff.png'))

%% Subfunctions
function [modelList, nModels] = getmodellist

    if isempty(getenv('IFILES'))
        error('IFILES environment variable not set')
    end

    modelSteffenFolder = fullfile(getenv('IFILES'), 'GIA', 'SteffenGrids');
    modelList = ls2cell(fullfile(modelSteffenFolder, 'Steffen*.mat'));
    % Remove the '_SD.mat' extension
    modelList = cellfun(@(x) x(1:end - 7), modelList, 'UniformOutput', false);

    nModels = length(modelList);
end

function GiaCorrAltimetry = getdisplacementdata(domainList, modelList, L)
    dataPath = mfilename('fullpath');
    dataPath = [dataPath, '-L', num2str(L), '.mat'];

    if exist(dataPath, 'file')
        warning('off', 'MATLAB:load:variableNotFound')
        load(dataPath, 'GiaCorrAltimetry')

        if exist('GiaCorrAltimetry', 'var')
            fprintf('%s loaded %s\n', upper(mfilename), dataPath)
            return
        end

    end

    nDomains = length(domainList);
    nModels = length(modelList);

    GiaCorrAltimetry = table('Size', [nModels, nDomains], ...
        'VariableTypes', [repmat({'double'}, [1, nDomains])], ...
        'VariableNames', domainList, ...
        'RowNames', modelList);

    parfor iModel = 1:nModels

        for iDomain = 1:nDomains
            domain = GeoDomain(domainList{iDomain}, "Buffer", 1, "DefaultParams", true); %#ok<PFBNS>
            [~, trend] = giaz2slept(modelList{iModel}, L, domain);
            trend = trend / domain.Area * 1e3;
            GiaCorrAltimetry{iModel, iDomain} = trend;
        end

    end

    try
        save(dataPath, 'GiaCorrAltimetry', '-append')
    catch
        save(dataPath, 'GiaCorrAltimetry')
    end

    fprintf('%s saved %s\n', upper(mfilename), dataPath)

end

function Altimetry = getaltimetrydata(domainList)
    dataPath = fullfile(getenv('IFILES'), 'OTHERS', 'Ludwigsen2024.xlsx');

    nDomains = length(domainList);
    domainSheetName = {'Global', 'North Pacific', 'South Pacific', 'North Atlantic', 'South Atlantic', 'Indian Ocean'};

    Altimetry = table('Size', [2, nDomains], ...
        'VariableTypes', [repmat({'double'}, [1, nDomains])], ...
        'VariableNames', domainList, 'RowNames', {'Alt-Str', 'GIA'});

    for iDomain = 1:nDomains
        AltimetryBasin = readtable(dataPath, "Sheet", domainSheetName{iDomain});
        f = polyfit(AltimetryBasin.time, AltimetryBasin.steAlt, 1);
        Altimetry{1, iDomain} = f(1);
        f = polyfit(AltimetryBasin.time, AltimetryBasin.GIA, 1);
        Altimetry{2, iDomain} = f(1);
        Altimetry{1, iDomain} = Altimetry{1, iDomain} + Altimetry{2, iDomain};
    end

end

function [Grace, GiaCorrMass] = getgracedata(domainList, modelList, L, timeRange)
    dataPath = mfilename('fullpath');
    dataPath = [dataPath, '-L', num2str(L), '.mat'];

    if exist(dataPath, 'file')
        warning('off', 'MATLAB:load:variableNotFound')
        load(dataPath, 'Grace', 'GiaCorrMass')

        if exist('Grace', 'var') && exist('GiaCorrMass', 'var')
            fprintf('%s loaded %s\n', upper(mfilename), dataPath)
            return
        end

    end

    nModels = length(modelList);
    nDomains = length(domainList);
    Grace = table('Size', [1, nDomains], ...
        'VariableTypes', [repmat({'double'}, [1, nDomains])], ...
        'VariableNames', domainList, ...
        'RowNames', {'MassTrend'});
    GiaCorrMass = table('Size', [nModels, nDomains], ...
        'VariableTypes', [repmat({'double'}, [1, nDomains])], ...
        'VariableNames', domainList, ...
        'RowNames', modelList);

    for iDomain = 1:nDomains
        domain = GeoDomain(domainList{iDomain}, 'Buffer', 1, 'DefaultParams', true);
        [date, slept] = grace2slept_new({'CSR', 'RL06', 60}, domain, 60, "Unit", 'SD', "TimeRange", timeRange, "BeQuiet", true);
        [~, ~, ~, ~, ~, ~, totalparams, ~, ~] = ...
            slept2resid_new(slept, date, [3, 365.0, 181.0], ...
            "Domain", domain, "Unit", 'year', "BeQuiet", true);

        Grace{1, iDomain} = mass2weq(totalparams(2), domain, 'seawater', 'mm');

        parfor iModel = 1:nModels
            model = modelList{iModel};

            if strcmp(model, 'ice6gd_vm5a')
                model = 'mascon';
            end

            [~, ~, ~, ~, total] = gia2slept(model, domain, L, "BeQuiet", true);
            % calSlept = slept - giaSlept;
            GiaCorrMass{iModel, iDomain} = mass2weq(total, domain, 'seawater', 'mm');
        end

    end

    try
        save(dataPath, 'Grace', 'GiaCorrMass', '-append')
    catch
        save(dataPath, 'Grace', 'GiaCorrMass')
    end

    fprintf('%s saved %s\n', upper(mfilename), dataPath)

end

function plottotalcomparison(Altimetry, Grace, modelList, domainDisplaynames)
    nModels = length(modelList);
    nDomains = length(domainDisplaynames);

    hold on

    for iModel = 1:nModels
        model = erase(modelList{iModel}, {'Steffen_', '-ice'});
        model = strsplit(model, '_');

        if contains(model{1}, 'ice6gd')
            mcolor = cc('g3');
            xOfs = -0.1;
        elseif contains(model{1}, 'ice6g')
            mcolor = cc('b3');
            xOfs = -0.1;
        elseif contains(model{1}, 'ice7g')
            mcolor = cc('r3');
            xOfs = 0;
        elseif contains(model{1}, 'anu')
            mcolor = cc('y3');
            xOfs = 0.1;
        else
            mcolor = 'y';
            xOfs = 0;
        end

        for iDomain = 1:nDomains
            xPlot = iDomain;

            scatter(xPlot - 0.2 + xOfs, ...
                Altimetry{iModel, iDomain}, ...
                'o', "MarkerEdgeColor", mcolor, ...
                "HandleVisibility", 'off')
            scatter(xPlot + 0.2 + xOfs, ...
                Grace{iModel, iDomain}, ...
                '+', "MarkerEdgeColor", mcolor, ...
                "HandleVisibility", 'off')
        end

    end

    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('g3'), "DisplayName", 'ice6g\_d')
    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('b3'), "DisplayName", 'ice6g\_c')
    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('r3'), "DisplayName", 'ice7g')
    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('y3'), "DisplayName", 'anu')
    scatter(nan, nan, "MarkerEdgeColor", 'none', "DisplayName", '')
    scatter(nan, nan, 'o', 'k', "DisplayName", 'Altimetry')
    scatter(nan, nan, '+', 'k', "DisplayName", 'GRACE')
    hold off

    legend("Box", 'off')
    ylabel('Ocean mass trend [mm/yr]')

    xticks(1:nDomains)
    xticklabels(domainDisplaynames)
    xlim([1, nDomains] + [-0.5, 0.5])
    box on
end

function plotgiacomparison(GiaCorr, modelList, domainDisplaynames, mshape)
    nModels = length(modelList);
    nDomains = length(domainDisplaynames);

    hold on

    for iModel = 1:nModels
        model = erase(modelList{iModel}, {'Steffen_', '-ice'});
        model = strsplit(model, '_');

        if contains(model{1}, 'ice6gd')
            mcolor = cc('g3');
            xOfs = -0.1;
        elseif contains(model{1}, 'ice6g')
            mcolor = cc('b3');
            xOfs = -0.1;
        elseif contains(model{1}, 'ice7g')
            mcolor = cc('r3');
            xOfs = 0;
        elseif contains(model{1}, 'anu')
            mcolor = cc('y3');
            xOfs = 0.1;
        else
            mcolor = 'y';
            xOfs = 0;
        end

        for iDomain = 1:nDomains
            xPlot = iDomain;

            scatter(xPlot + xOfs, ...
                GiaCorr{iModel, iDomain}, ...
                mshape, "MarkerEdgeColor", mcolor, ...
                "HandleVisibility", 'off')
        end

    end

    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('g3'), "DisplayName", 'ice6g\_d')
    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('b3'), "DisplayName", 'ice6g\_c')
    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('r3'), "DisplayName", 'ice7g')
    scatter(nan, nan, 'o', "MarkerEdgeColor", cc('y3'), "DisplayName", 'anu')
    scatter(nan, nan, "MarkerEdgeColor", 'none', "DisplayName", '')
    scatter(nan, nan, 'o', 'k', "DisplayName", 'Altimetry')
    scatter(nan, nan, '+', 'k', "DisplayName", 'GRACE')
    hold off

    yline(0, 'k--', "HandleVisibility", 'off')
    ylim([-6, 1])

    legend("Box", 'off', "Location", 'southeast')
    ylabel('GIA trend [mm/yr]')

    xticks(1:nDomains)
    xticklabels(domainDisplaynames)
    xlim([1, nDomains] + [-0.5, 0.5])
    box on
end
