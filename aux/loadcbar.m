%% LOADCBAR
% Loads a colourbar with specified settings.
% This function loads a colourbar with specified settings. It creates a
% colourmap using the specified colourmap name and number of colour levels.
% The colourbar is then adjusted to the specified colour limits and
% labelled with the specified name. Additional figure settings are provided
% through the FSettings structure.
%
% Syntax
%   loadcbar(levels)
%   loadcbar(limits, h)
%   loadcbar(__, 'Title', title, 'CmapName', name)
%   loadcbar(__, 'FSettings', FigureSettings)
%   [cbar, levels] = loadcbar(__)
%
% Input arguments
%  levels - Colour levels
%  limits - Step size for colour levels
%  h - Step size for colour levels
%  title - Name of the colourbar
%  name - Name of the colourmap
%  FigureSettings - Structure with the figure settings
%      If given, it should have the following fields:
%      - FsizeCaption: Font size for the caption
%      - Fname: Font name
%      - LwidthThin: Line width for thin lines
%      The default value is [], which means that the default MATLAB
%      settings are used.
%
% Output arguments
%   - cbar - Handle to the colourbar
%   - cLevels - Colour levels
%
% Example
%   >>  cbarLabel = 'Temperature anmomaly [°C]';
%   >>  cStep = 0.5;
%   >>  cLim = [0, 10];
%   >>  cmapName = 'inverted temperature anomaly';
%   >>  FSettings.LwidthThick = 2;
%   >>  FSettings.Fname = 'Arial';
%   >>  FSettings.FsizeCaption = 12;
%   >>  [cLevels, cbar] = loadcbar(cLim, cStep, ...
%       'Title', cbarLabel, 'CmapName', cmapName, FSettings)
%
% Last modified by
%   2021/08/13, williameclee@arizona.edu (@williameclee)
%   2021/07/02, williameclee@arizona.edu (@williameclee)

function varargout = loadcbar(varargin)
    %% Initialisation
    % Parse inputs
    [cLevels, cLim, cbarLabel, cmapName, FSettings, location] = ...
        parseinputs(varargin{:});

    %% Compute the colour levels and centre
    [~, cLevels] = loadcmap(cmapName, cLevels);
    cLevelsForPlot = coarseclevels(cLevels);

    %% Create and style the colourbar
    cbar = colorbar(location);
    cbar.Ticks = cLevelsForPlot;
    clim(cLim)

    if ~isempty(cbarLabel)
        cbar.Label.String = cbarLabel;
    end

    if ~isempty(FSettings)
        set(cbar, 'LineWidth', FSettings.LwidthThin, ...
            'FontName', FSettings.Fname, ...
            'FontSize', FSettings.FsizeNormal, 'FontWeight', 'normal')
    end

    if nargout == 0
        return
    end

    varargout = {cbar, cLevels};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addOptional(p, 'cLim', [-1, 1], @isnumeric);
    addOptional(p, 'cStep', [], @isnumeric);
    addOptional(p, 'Title', '', @ischar);
    addOptional(p, 'Colormap', 'seismic', @(x) ischar(x) || isstring(x));
    addOptional(p, 'FSettings', [], @isstruct);
    addParameter(p, 'Location', 'eastoutside', @ischar);
    parse(p, varargin{:});

    cLim = p.Results.cLim;
    cStep = p.Results.cStep;

    if isscalar(cLim)
        error('cLim must be a 2-element numeric array.');
    elseif length(cLim) == 2

        if isempty(cStep)
            cStep = (cLim(2) - cLim(1)) / 16;
        end

        cLevels = cLim(1):cStep:cLim(2);
    else
        cLevels = cLim;
    end

    cLim = [min(cLevels), max(cLevels)];
    cbarTitle = p.Results.Title;
    cmapName = char(p.Results.Colormap);
    FSettings = p.Results.FSettings;
    location = p.Results.Location;

    varargout = {cLevels, cLim, cbarTitle, cmapName, FSettings, location};
end

function cLevelsC = coarseclevels(cLevels)

    if length(cLevels) < 10
        cLevelsC = cLevels;
        return
    end

    zeroId = find(cLevels == 0);

    nSkip = floor(length(cLevels) / 9) + 1;

    if nSkip == 1
        cLevelsC = cLevels;
        return
    end

    if isempty(zeroId) || nmod(zeroId, nSkip) == 1
        cLevelsC = cLevels(1:nSkip:end);
    else
        cLevelsC = cLevels(nmod(zeroId, nSkip):nSkip:zeroId - 1);
        cLevelsC = [cLevelsC, cLevels(zeroId:nSkip:end)];
    end

end

function x = nmod(x, n)
    x = mod(x, n);
    x(x == 0) = n;
end
