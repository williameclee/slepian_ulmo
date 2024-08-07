function dataFileAttr = dataattrchar(varargin)
    p = inputParser;
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 0);
    addOptional(p, 'Latlim', 90, ...
        @(x) isnumeric(x) && length(x) <= 2);
    addOptional(p, 'MoreBuffers', []);
    parse(p, varargin{:});
    upscale = p.Results.Upscale;
    latlim = p.Results.Latlim;
    buf = p.Results.Buffer;
    moreBuf = p.Results.MoreBuffers;

    if length(latlim) == 2

        if latlim(1) == -latlim(2)
            latlim = max(latlim);
        end

    end

    %% Find the data file
    dataFileAttr = cell(1, 3);

    if ~isempty(upscale)

        if ~(upscale == 0 || upscale == 1)
            dataFileAttr{1} = num2str(upscale);
        end

    end

    if ~(buf == 0) || ~isempty(moreBuf)
        dataFileAttr{2} = num2str(buf);

        if ~isempty(moreBuf)
            moreBuf = formatmorebuf(moreBuf);

            for i = 1:2:length(moreBuf)
                moreBuf{i} = [moreBuf{i}, num2str(moreBuf{i + 1})];
            end

            moreBuf = char(join(moreBuf(1:2:end), '_'));

            dataFileAttr{2} = [dataFileAttr{2}, '_', moreBuf];
        end

    end

    if any(isnan(latlim)) || any(isempty(latlim))
    elseif isscalar(latlim)

        if ~(latlim == 90 || isempty(latlim))
            dataFileAttr{3} = num2str(latlim);
        end

    elseif ~isequal(latlim, [-90, 90])
        inclangSign = [];

        if latlim(1) < 0
            inclangSign(1) = 's';
        elseif latlim(1) > 0
            inclangSign(1) = 'n';
        end

        if latlim(2) < 0
            inclangSign(2) = 's';
        elseif latlim(2) > 0
            inclangSign(2) = 'n';
        end

        dataFileAttr{3} = ...
            [num2str(abs(latlim(1))), inclangSign(1), '_', ...
             num2str(abs(latlim(2))), inclangSign(2)];

    end

    % Find the last filled attribute
    iLastUsedAttr = find(~cellfun(@isempty, dataFileAttr), 1, 'last');

    if ~isempty(iLastUsedAttr)
        dataFileAttr = dataFileAttr(1:iLastUsedAttr);
        % Fill in 0
        emptyBeforeLast = 1:(iLastUsedAttr - 1);
        emptyCellsIndex = cellfun(@isempty, dataFileAttr(emptyBeforeLast));
        dataFileAttr(emptyBeforeLast(emptyCellsIndex)) = {'0'};
        % Add a hyphen before each attribute
        if length(dataFileAttr) > 1
            dataFileAttr(2:end) = cellfun(@(x) ['-', x], dataFileAttr(2:end), ...
                'UniformOutput', false);
        end

    end

    if all(cellfun(@isempty, dataFileAttr))
        dataFileAttr = '';
    else
        dataFileAttr = char(join(dataFileAttr, ''));
    end

end

%% Subfunctions
function moreBuf = formatmorebuf(moreBuf)
    moreBufDomain = lower(moreBuf(1:2:end));
    moreBufWidth = moreBuf(2:2:end);
    [moreBuf(1:2:end), moreBufSortId] = sort(moreBufDomain);
    moreBuf(2:2:end) = moreBufWidth(moreBufSortId);

    for i = 1:2:length(moreBuf)

        if strcmp(moreBuf{i}, 'earthquakes')
            moreBuf{i} = 'eqs';
        elseif strcmp(moreBuf{i}, 'africa')
            moreBuf{i} = 'AF';
        elseif strcmp(moreBuf{i}, 'antarctica')
            moreBuf{i} = 'AN';
        elseif strcmp(moreBuf{i}, 'australia')
            moreBuf{i} = 'AU';
        elseif strcmp(moreBuf{i}, 'eurasia')
            moreBuf{i} = 'EA';
        elseif strcmp(moreBuf{i}, 'greenland')
            moreBuf{i} = 'GL';
        elseif strcmp(moreBuf{i}, 'namerica')
            moreBuf{i} = 'NA';
        elseif strcmp(moreBuf{i}, 'samerica')
            moreBuf{i} = 'SA';
        elseif strcmp(moreBuf{i}, 'npacific')
            moreBuf{i} = 'npac';
        elseif strcmp(moreBuf{i}, 'spacific')
            moreBuf{i} = 'spac';
        elseif strcmp(moreBuf{i}, 'natlantic')
            moreBuf{i} = 'natl';
        elseif strcmp(moreBuf{i}, 'satlantic')
            moreBuf{i} = 'satl';
        end

        moreBuf{i} = capitalise(moreBuf{i});

    end

end
