%% LOADCMAP
% Loads a colourmap with specified name and colour levels.
%
% Syntax
%    loadcmap
%    loadcmap(Name)
%    loadcmap(Name, levels)
%    loadcmap(Name, limits, h)
%    [cmap, levels, limits] = loadcmap(__)
%
% Input arguments
%   Name - Name of the colourmap
%        The default value is 'seismic', which is a custom colourmap (see notes below).
%   levels - Levels of the colourmap
%   limits - Limits of the colourmap
%   h - Spacing between the levels
%
% Output arguments
%   cmap - N-by-3 colourmap matrix of RGB values
%   levels - Levels of the colourmap
%   limits - Limits of the colourmap
%
% Notes
%   The function tries to call custom colourmaps from the 'MatlabColourmapGenerator' repository:
%   https://github.com/williameclee/MatlabColourmapGenerator.git
%   If the colourmap is not available, it will default to the 'jet' colourmap.
%
% Examples
%   Load the seismic colourmap with default levels and limits
%   >> loadcmap
%   Load the seismic colourmap with custom levels
%   >> loadcmap('seismic', linspace(-1, 1, 32))
%
% See also
%   LOADCBAR, KEYNOTECMAP (KCMAP), JET, COLORMAP
%
% Last modified by
%   2021/08/13, williameclee@arizona.edu (@williameclee)

function varargout = loadcmap(varargin)
    p = inputParser;
    addOptional(p, 'Name', 'seismic', @ischar);
    addOptional(p, 'cLim', [-1, 1], @isnumeric);
    addOptional(p, 'cStep', 16, @isnumeric);
    parse(p, varargin{:});
    cmapName = p.Results.Name;
    cLim = p.Results.cLim;
    cStep = p.Results.cStep;

    %% Compute the colour levels and centre
    if isscalar(cLim)
        error('Invalid clim input.');
    elseif length(cLim) == 2
        cLevels = cLim(1):cStep:cLim(2);
    else
        cLevels = cLim;
    end

    %% Load the colourmap
    try
        cmap = ccmap(cmapName, cLevels);
    catch

        try
            cmap = kcmap(cmapName, cLevels);
        catch
            cmap = jet(length(cLevels));
            fprintf('The colourmap %s is not available. \n', cmapName);
            fprintf('Consider getting better colourmaps at: \n%s\n', ...
            'https://github.com/williameclee/MatlabColourmapGenerator.git');
        end

    end

    cmap = max(cmap, 1e-3);
    colormap(cmap)

    try
        clim(cLim)
    catch
    end

    if nargout == 0
        clear varargout
        return
    end

    cLim = [min(cLevels), max(cLevels)];
    varargout = {cmap, cLevels, cLim};

end
