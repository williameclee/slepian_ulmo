%% EXAMPLE: GIA COMPARISON
% Compares how different GIA models affect the trend in mass in different ocean basins.

function giacomparison_example
    L = 60;
    timeRange = datetime([2003, 2022], [1, 12], [1, 31]);

    domainList = {'oceans', 'npacific', 'spacific', 'natlantic', 'satlantic', 'indian'};
    nDomains = length(domainList);

    [modelList, nModels] = getmodellist;

    trend = zeros([nDomains, nModels]);

    for iDomain = 1:nDomains
        domain = GeoDomain(domainList{iDomain}, 'Buffer', 1, 'DefaultParams', true);

        [date, slept] = grace2slept_new({'CSR', 'RL06', 60}, domain, 60, "Unit", 'SD', "TimeRange", timeRange, "BeQuiet", true);

        parfor iModel = 1:nModels
            model = modelList{iModel};
            [~, giaSlept] = gia2slept(date, model, domain, L, "BeQuiet", true);
            sleptCor = slept - giaSlept;
            [~, ~, ~, ~, ~, ~, totalparams, ~, ~] = ...
                slept2resid_new(sleptCor, date, [3, 365.0, 181.0], ...
                "Domain", domain, "Unit", 'year', "BeQuiet", true);

            trend(iDomain, iModel) = mass2weq(totalparams(2), domain, 'seawater', 'mm');

            if iscell(model)
                modelStr = strjoin(model, '_');
            else
                modelStr = model;
            end

            fprintf('Trend in %s is %.2f mm/year (GIA: %s)\n', ...
                domain.DisplayName, trend(iDomain, iModel), modelStr);
        end

    end

    plotgiacomparison(timeRange, domainList, trend, modelList)
end

%% Subfunctions
function [modelList, nModels] = getmodellist
    modelSteffenFolder = fullfile(getenv('IFILES'), 'GIA', 'SteffenGrids');
    modelSteffen = ls2cell(fullfile(modelSteffenFolder, 'Steffen*.mat'));
    % Remove the '_SD.mat' extension
    modelSteffen = cellfun(@(x) x(1:end - 7), modelSteffen, 'UniformOutput', false);
    modelList = ['Paulson07', modelSteffen, 'mascon'];

    nModels = length(modelList);
end

function plotgiacomparison(timeRange, domainList, trend, modelList)
    modelGroups = {'nu-ice', 'ice6g', 'ice7g', 'mascon', 'Paulson07'};
    modelColours = {kc('b3'), kc('g3'), kc('c3'), kc('y3'), kc('r3')};
    modelXoffset = [-0.1, 0.1, 0.2, 0, 0];

    figure(1)
    clf

    for iDomain = 1:length(domainList)
        hold on
        boxchart(ones(size(trend(iDomain, :))) * iDomain, trend(iDomain, :), ...
            "BoxFaceColor", kc('k3'), "BoxEdgeColor", kc('k3'), "WhiskerLineColor", kc('k3'), ...
            "HandleVisibility", 'off');
        hold off

        for iGroup = 1:length(modelGroups)
            trend2plot = trend(iDomain, contains(modelList, modelGroups{iGroup}));
            x2plot = iDomain + modelXoffset(iGroup);
            hold on
            scatter(x2plot, trend2plot, ...
                'MarkerFaceColor', modelColours{iGroup}, 'MarkerEdgeColor', modelColours{iGroup}, ...
                "HandleVisibility", 'off');
            hold off
        end

    end

    for iGroup = 1:length(modelGroups)
        hold on
        scatter(nan, nan, ...
            'MarkerFaceColor', modelColours{iGroup}, 'MarkerEdgeColor', modelColours{iGroup}, ...
            "DisplayName", modelGroups{iGroup});
        hold off
    end

    legend("Location", 'southwest', "Box", 'off')

    xlim([1, length(domainList)] + [-1, 1] * 0.5)
    xticks(1:length(domainList))
    xticklabels(cellfun(@(x) GeoDomain(x).DisplayName, domainList, 'UniformOutput', false))

    xlabel('Ocean basins')
    ylabel('Trend [mm/year]')
    title(sprintf('Mass in ocean basins with different GIA models, %i-%i', ...
        year(timeRange(1)), year(timeRange(2))))

    set(gca, "Box", 'on', "Layer", 'top')
    set(gcf, "NumberTitle", 'off', "Name", ...
        sprintf('Mass in ocean basins with different GIA models (%s)', upper(mfilename)))

end
