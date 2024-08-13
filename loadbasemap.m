%% LOADBASEMAP
% Loads a axesm-based axes object for a given geographic domain
%
% Syntax
%   loadbasemap(domain)
%   loadbasemap(domain, FigureSettings)
%   loadbasemap(domain, FigureSettings, latlim, lonlim)
%   loadbasemap(__, 'PrjMethod', prjMethod)
%   loadbasemap(__, 'MlabelPos', mlabelPos)
%   [axm, latlim, lonlim] = loadbasemap(__)
%
% Input arguments
%   domain - Geographic domain or a latitude-longitude pair
%       - A geographic domain (GeoDomain object).
%       - A string of the domain name.
%       - A cell array of the form {'domain name', buf}.
%   FigureSettings - Structure with the figure settings
%       If given, it should have the following fields:
%       - FsizeCaption: Font size for the caption
%       - Fname: Font name
%       - LwidthThin: Line width for thin lines
%       The default value is [], which means that the default MATLAB 
%       settings are used.
%   latlim - Latitude limits
%   lonlim - Longitude limits
%   prjMethod - Projection method
%       It should preferably be area-preserving.
%       The default value is 'lambcyln' for most regions, and 'wagner4' 
%       for oceans.
%   mlabelPos - Meridian label position
%       The default value is 'south' for most regions, and 'north' for 
%       Antarctica.
%
% Output arguments
%   axm - Axes object
%   latlim - Latitude limits
%   lonlim - Longitude limits
%
% Examples
%   Load the basemap for the world
%   >> loadbasemap
%   Load the basemap for the world with a custom figure settings
%   >> loadbasemap('world', struct('LwidthThin', 0.5, ...
%       'FsizeCaption', 12, 'Fname', 'Arial'))
%   Load the basemap for the world with custom limits
%
% Last modified by
%   2024/08/13, williameclee@arizona.edu (@williameclee)

function varargout = loadbasemap(varargin)
    %% Initialisation
    % Get manually defined values
    [domain, FSettings, latLimM, lonLimM, prjMethodM, mlabelPosM] = ...
        parseinputs(varargin{:});
    % Fallback values
    prjMethod = 'lambcyln';
    mlabelPos = 'south';
    % Region-specific values
    [latLim, lonLim] = loaddefaultregionlimits(domain);

    switch domain
        case 'antarctica'
            prjMethod = 'eqaazim';
            mlabelPos = 'north';
        case 'greenland'
            prjMethod = 'eqaconic';
            mlabelPos = 'south';
    end

    if isocean(domain)
        prjMethod = 'wagner4';
    end

    % if strcmp(region, 'oceans')
    %     lonLim = [20, 380];
    % end

    latLim = conddefval(latLimM, latLim);
    lonLim = conddefval(lonLimM, lonLim);
    prjMethod = conddefval(prjMethodM, prjMethod);
    mlabelPos = conddefval(mlabelPosM, mlabelPos);

    AxesmDefaultSettings = ...
        {'Frame', 'on', 'Grid', 'off', ...
         'MeridianLabel', 'off', 'ParallelLabel', 'off', ...
         'MLabelParallel', mlabelPos, ...
         'GLineStyle', ':'};

    switch prjMethod
        case 'eqaazim'
            axm = axesm('MapProjection', prjMethod, ...
                'Origin', [mean(latLim), mean(lonLim), 0], ...
                AxesmDefaultSettings{:});
        otherwise
            axm = axesm('MapLatLimit', latLim, 'MapLonLimit', lonLim, ...
                'MapProjection', prjMethod, ... % area-preserving
                AxesmDefaultSettings{:});
    end

    if ~isempty(FSettings)
        set(axm, 'LineWidth', FSettings.LwidthThin * 2, ...
            'FontSize', FSettings.FsizeCaption, ...
            'FontName', FSettings.Fname);
        set(axm, 'GridLineWidth', FSettings.LwidthThin);
    end

    try
        set(axm, 'GridColor', kc('k3'));
    catch
        set(axm, 'GridColor', 'k');
    end

    tightmap
    set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')

    if nargout == 0
        return
    end

    varargout = {axm, latLim, lonLim};
end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addOptional(p, 'Domain', 'world', ...
        @(x) ischar(x) || isstring(x) || iscell(x) || isa(x, "GeoDomain"));
    addOptional(p, 'FigureSettings', [], ...
        @isstruct);
    addOptional(p, 'LatLim', [], ...
        @(x) (isnumeric(x) && length(x) == 2) || isempty(x));
    addOptional(p, 'LonLim', [], ...
        @(x) (isnumeric(x) && length(x) == 2) || isempty(x));
    addParameter(p, 'PrjMethod', [], ...
        @(x) ischar(x) || isstring(x) || isempty(x));
    addParameter(p, 'MlabelPos', [], ...
        @(x) ischar(x) || isnumeric(x) || isempty(x));
    parse(p, varargin{:});

    domain = p.Results.Domain;

    if iscell(domain)
        domain = domain{1};
    elseif isa(domain, "GeoDomain")
        domain = domain.Domain;
    end

    FSettings = p.Results.FigureSettings;
    latLim = p.Results.LatLim;
    lonLim = p.Results.LonLim;
    prjMethod = p.Results.PrjMethod;
    mlabelPos = p.Results.MlabelPos;

    varargout = {domain, FSettings, latLim, lonLim, prjMethod, mlabelPos};
end
