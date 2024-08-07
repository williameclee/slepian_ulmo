%% SPACIFIC
% Finds the longitude and latitude coordinates of the South Pacific Ocean.
%
% Syntax
%   XY = spacific(upscale, buf)
%   XY = spacific(__, latlim)
%   XY = spacific(__, morebuffers)
%   XY = spacific(__, 'Name', value)
%   spacific(__)
%
% Input arguments
%   upscale - How many times to upscale the data
%       This option should be used carefully, as spline interpolation do
%       not work well at sharp corners, and the original sampling rate
%       should be sufficient for most applications.
%       The default value is 0 (no upscaling).
%   buffer - The (negative) buffer from the coastlines in degrees
%       The default value is 0 (no buffer).
%   latlim - The latitudes of the polar caps in degrees
%       The inclination angle can be specified as a scalar or a 1-by-2
%       vector.
%       - If a scalar, the inclination angle is applied to both the north
%           and south polar caps.
%       - If a 1-by-2 vector, the first element is the inclination angle of
%           the north polar cap, and the second element is the inclination
%           angle of the south polar cap.
%       The default value is 90 (no polar caps).
%   morebuffers - Additional buffers to apply to the coastlines
%       The additional buffers must be specified as a cell array of
%       domain names and buffer values.
%       The default value is an empty cell array.
%   LonOrigin - The longitude origin of the data
%       The domain will be contained within the range
%       [LonOrigin - 180, LonOrigin + 180].
%       The default value is 180 (i.e. the range of longitude is [0, 360]).
%   ForceNew - Force the function to reload the data
%       The default value is false.
%   BeQuiet - Suppress the output messages
%       The default value is false.
%
% Output arguments
%   XY - Closed-curved coordinates of the South Pacific Ocean
%
% Examples
%   The 'default' ocean domain used for most applications is given by
%   XY = spacific('Buffer', 1, 'Latlim', 60);
%
% Data source
%   The coastline data is based on the GSHHG coastline data:
%       Wessel, P., & Smith, W. H. F. (1996).
%       doi: 10.1029/96JB00104
%   The ocean boundaries are based on IHO's 'Limits of oceans and seas':
%       International Hydrographic Organization & Sieger, R. (2012).
%       doi: 10.1594/PANGAEA.777975
%
% Last modified by
%   williameclee-at-arizona.edu, 2024/08/07

function varargout = spacific(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, latlim, buf, moreBufs, lonOrigin, forceNew, beQuiet] = ...
        parsecoastinputs(varargin);
    oceanParts = ...
        {'South Pacific Ocean, eastern part', ...
     'South Pacific Ocean, western part'};
    % Title for plots
    oceanTitle = sprintf('%s (%s)', domainname(mfilename, 'long'), upper(mfilename));

    %% Check if the data file exists
    [dataFile, ~, dataFileExists] = coastfilename(mfilename, ...
        'Upscale', upscale, 'Latlim', latlim, ...
        'Buffer', buf, 'MoreBuffers', moreBufs);

    if dataFileExists && ~forceNew
        load(dataFile, 'XY')

        if beQuiet < 2
            fprintf('%s loaded %s\n', upper(mfilename), dataFile)
        end

        varargout = returncoastoutputs(nargout, XY, oceanTitle);

        return
    end

    %% Compute the data
    [oceanPoly, oceanLatlim, oceanLonlim] = ...
        findoceanboundary(oceanParts, latlim, lonOrigin);

    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBufs, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    oceanPoly = subtract(oceanPoly, coastPoly);

    XY = closecoastline(oceanPoly.Vertices);

    % Convert the data to the right format
    [X, Y] = poly2cw(XY(:, 1), XY(:, 2));
    XY = removeduplicatevertices([X, Y]);
    % Make sure the minimum longitude is greater than 0
    XY(:, 1) = XY(:, 1) - floor(min(XY(:, 1)) / 360) * 360;

    %% Save and return required data
    save(dataFile, 'XY', '-v7.3')

    if beQuiet < 2
        fprintf('%s saved %s\n', upper(mfilename), dataFile)
    end

    varargout = returncoastoutputs(nargout, XY, oceanTitle);

end

function coastPoly = manualadjustment(coastPoly, buf, lonOrigin)

    coastPoly = addlandregion(coastPoly, ...
        'Latlim', [-53, -50; -54, -53], ...
        'Lonlim', [289, 292; 291, 292], ...
        'LongitudeOrigin', lonOrigin);
    % coastPoly = addlandregion(coastPoly, ...
    %     'Latlim', [-54, -53], 'Lonlim', [291, 292], ...
    %     'LongitudeOrigin', lonOrigin);

    if buf >= 0.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-1, 0; -11, -9; -41, -40; -41, -39; -77, -76], ...
            'Lonlim', [130, 132; 141, 142.5; 172.5, 175; 173.5, 174; 164, 166], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-40, -39; -3, 0], 'Lonlim', [143.5, 145.5; 134, 138], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 3.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-3, 0], 'Lonlim', [144, 152], ...
            'LongitudeOrigin', lonOrigin);
    end

end
