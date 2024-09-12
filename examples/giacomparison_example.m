%% EXAMPLE: GIA COMPARISON
% Compares how different GIA models affect the trend in mass in different ocean basins.

function [trend, giaTrend] = giacomparison_example
    L = 60;
    timeRange = datetime([2003, 2022], [1, 12], [1, 31]);

    domainList = {'oceans', 'npacific', 'spacific', 'natlantic', 'satlantic', 'indian'};
    nDomains = length(domainList);

    [modelList, nModels] = getmodellist;

    dataPath = mfilename('fullpath');
    dataPath = [dataPath, '.mat'];

    if exist(dataPath, 'file')
        load(dataPath, 'trend', 'giaTrend')
    else
        trend = zeros([nDomains, nModels]);
        giaTrend = zeros([nDomains, nModels]);

        for iDomain = 1:nDomains
            domain = GeoDomain(domainList{iDomain}, 'Buffer', 1, 'DefaultParams', true);

            [date, slept] = grace2slept_new({'CSR', 'RL06', 60}, domain, 60, "Unit", 'SD', "TimeRange", timeRange, "BeQuiet", true);

            parfor iModel = 1:nModels
                model = modelList{iModel};
                [~, giaSlept] = gia2slept(date, model, domain, L, "BeQuiet", true);
                calSlept = slept - giaSlept;
                [~, ~, ~, ~, ~, ~, totalparams, ~, ~] = ...
                    slept2resid_new(calSlept, date, [3, 365.0, 181.0], ...
                    "Domain", domain, "Unit", 'year', "BeQuiet", true);

                trend(iDomain, iModel) = mass2weq(totalparams(2), domain, 'seawater', 'mm');
                [~, ~, ~, ~, ~, ~, totalparams, ~, ~] = ...
                    slept2resid_new(giaSlept, date, [3, 365.0, 181.0], ...
                    "Domain", domain, "Unit", 'year', "BeQuiet", true);

                giaTrend(iDomain, iModel) = mass2weq(totalparams(2), domain, 'seawater', 'mm');

                if iscell(model)
                    modelStr = strjoin(model, '_');
                else
                    modelStr = model;
                end

                fprintf('Total trend in %s is %.2f mm/year (GIA model: %s)\n', ...
                    domain.DisplayName, trend(iDomain, iModel), modelStr);
                fprintf('GIA trend in %s is %.2f mm/year (GIA model: %s)\n', ...
                    domain.DisplayName, trend(iDomain, iModel), modelStr);
            end

        end

        save(dataPath, 'trend', 'giaTrend')

    end

    plotgiacomparison(timeRange, domainList, trend, modelList)
end

%% Subfunctions
function [modelList, nModels] = getmodellist

    if isempty(getenv('IFILES'))
        error('IFILES environment variable not set')
    end

    modelSteffenFolder = fullfile(getenv('IFILES'), 'GIA', 'SteffenGrids');
    modelSteffen = ls2cell(fullfile(modelSteffenFolder, 'Steffen*.mat'));
    % Remove the '_SD.mat' extension
    modelSteffen = cellfun(@(x) x(1:end - 7), modelSteffen, 'UniformOutput', false);
    modelList = ['Paulson07', modelSteffen, 'mascon'];

    nModels = length(modelList);
end

function plotgiacomparison(~, domainList, trend, modelList)
    modelGroups = {'anu-ice', 'ice6g', 'ice7g', 'mascon', 'Paulson07'};
    modelColours = {'b', 'g', 'c', 'y', 'r'};
    modelXoffset = [0, 0.1, 0.2, -0.2, -0.2];

    xlimRange = [1, length(domainList)] + [-1, 1] * 0.5;

    figure(1)
    clf

    hold on
    plot(xlimRange, [0, 0], ...
        "Color", kc('w3'), "LineStyle", '--', "HandleVisibility", 'off')
    hold off

    for iDomain = 1:length(domainList)
        hold on
        boxchart(ones(size(trend(iDomain, :))) * iDomain, trend(iDomain, :), ...
            "BoxFaceColor", kc('k1'), "BoxEdgeColor", kc('k2'), "WhiskerLineColor", kc('k2'), ...
            "HandleVisibility", 'off');
        hold off

        for iModel = 1:length(modelList)
            model = modelList{iModel};

            % if contains(model, 'ice7g')
            %     continue
            % end

            msize = 60;
            shape = 'o';
            shade = 3;

            if contains(model, 'Steffen')
                modelspec = strsplit(model, '_');

                if length(modelspec{3}) == 3 && ~contains(modelspec{3}, 'vm')

                    if strcmp(modelspec{3}(end), '2')
                        shape = 'diamond';
                    elseif strcmp(modelspec{3}(end), '6')
                        shape = 'square';
                    else
                        error('Unknown Steffen model: %s', modelspec{3})
                    end

                    switch modelspec{3}(1)
                        case 'f'
                            shade = 1;
                        case 'i'
                            shade = 2;
                        case 'l'
                            shade = 3;
                        case 'o'
                            shade = 4;
                    end

                end

            end

            groupId = find(cellfun(@(x) contains(model, x), modelGroups));
            trend2plot = trend(iDomain, iModel);
            x2plot = iDomain + modelXoffset(groupId);

            % if strcmp(shape, 'diamond')
            %     x2plot = x2plot + 0.02;
            % elseif strcmp(shape, 'square')
            %     x2plot = x2plot - 0.02;
            if contains(model, 'vm5a')
                x2plot = x2plot -0.3;
            elseif contains(model, 'vm7')
                x2plot = x2plot -0.4;
            end

            color = modelColours{groupId};
            hold on
            scatter(x2plot, trend2plot, shape, "SizeData", msize, ...
                'MarkerFaceColor', kc(color, shade), 'MarkerEdgeColor', 'k', ...
                "HandleVisibility", 'off');
            hold off
        end

    end

    hold on

    for iGroup = 1:length(modelGroups)

        if contains(modelGroups{iGroup}, 'ice7g')
            continue
        end

        scatter(nan, nan, 'o', "SizeData", msize, ...
            'MarkerFaceColor', kc(modelColours{iGroup}, 2), 'MarkerEdgeColor', 'k', ...
            "DisplayName", modelGroups{iGroup});
    end

    for i = 1:2
        scatter(nan, nan, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w', ...
            "DisplayName", '');
    end

    scatter(nan, nan, 'diamond', "SizeData", msize, ...
        'MarkerFaceColor', kc('w', 2), 'MarkerEdgeColor', 'k', ...
        "DisplayName", 'Low viscosity');
    scatter(nan, nan, 'square', "SizeData", msize, ...
        'MarkerFaceColor', kc('w', 2), 'MarkerEdgeColor', 'k', ...
        "DisplayName", 'High viscosity');

    lithThickness = [60, 90, 120, 150];

    for iShade = 1:4
        scatter(nan, nan, 'o', "SizeData", msize, ...
            'MarkerFaceColor', kc('w', iShade), 'MarkerEdgeColor', 'k', ...
            "DisplayName", sprintf('%i km lithosphere', lithThickness(iShade)));
    end

    hold off

    legend("Location", 'southwest', "Box", 'off', "NumColumns", 2);

    xlim(xlimRange)
    xticks(1:length(domainList))
    xticklabels(cellfun(@(x) GeoDomain(x).DisplayName, domainList, 'UniformOutput', false))

    xlabel('Ocean basins')
    ylabel('Trend [mm/year]')
    % title(sprintf('Mass in ocean basins by different GIA models, %i-%i', ...
    %     year(timeRange(1)), year(timeRange(2))))

    set(gca, "Box", 'on', "Layer", 'top')
    set(gcf, "NumberTitle", 'off', "Name", ...
        sprintf('Mass in ocean basins by different GIA models (%s)', upper(mfilename)))

    % %% Style
    % set(gca, "Box", 'off', "Layer", 'bottom')
    % set(gca, "FontName", 'Fira Code', "FontSize", 10)
    % set(gca, "Color", [26, 30, 75] / 255)
    % set(gca, "LineWidth", 1)
    % set(gcf, "Color", [0, 0, 0])
    % legend("TextColor", 'w', "Color", [26, 30, 75] / 255, "Box", 'on', "EdgeColor", [26, 30, 75] / 255)
    % set(gca, "XColor", [102, 153, 255] / 255, "YColor", [102, 153, 255] / 255, "GridColor", [102, 153, 255] / 255)
    % ax = gca;
    % ax.XLabel.Color = 'w';
    % ax.YLabel.Color = 'w';
    % ax.TitleFontSizeMultiplier = 1;
    % set(gca, "YGrid", 'on', "GridLineWidth", 1, "GridColor", 'k', "GridAlpha", 1)

    % for i = 1:length(ax.XTickLabel)
    %     ax.XTickLabel{i} = sprintf('\\color{white}%s', ax.XTickLabel{i});
    % end

    % for i = 1:length(ax.YTickLabel)
    %     ax.YTickLabel{i} = sprintf('\\color{white}%s', ax.YTickLabel{i});
    % end

    % annotation('rectangle', [0, 0, 1, 1], "Color", 'k', "LineWidth", 0.05)
    % set(gca, "Position", [120, 120, 960 - 180, 960 - 180] / 960)
    % set(gcf, "Units", 'pixels', "Position", [100, 100, 960, 960] / 2)

    % % set(gcf, "PaperUnits", 'inches', "PaperSize", [960, 960] / 72)

    % exportgraphics(gcf, [mfilename, '.pdf'], "ContentType", 'vector', "BackgroundColor", 'none')

end
