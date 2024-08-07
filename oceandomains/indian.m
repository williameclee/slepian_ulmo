%% INDIAN
% Finds the longitude and latitude coordinates of the Indian Ocean.
%
% Syntax
%   XY = indian(upscale, buf)
%   XY = indian(__, latlim)
%   XY = indian(__, morebuffers)
%   XY = indian(__, 'Name', value)
%   indian(__)
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
%   XY - Closed-curved coordinates of the Indian Ocean
%
% Examples
%   The 'default' ocean domain used for most applications is given by
%   XY = indian('Buffer', 1, 'Latlim', 60, ...
%       'MoreBuffers', {'earthquakes', 10});
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

function varargout = indian(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, latlim, buf, moreBuf, lonOrigin, forceNew, beQuiet] = ...
        parsecoastinputs(varargin);
    oceanParts = 'Indian Ocean';
    % Title for plots
    oceanTitle = sprintf('%s (%s)', domainname(mfilename, 'long'), upper(mfilename));

    %% Check if the data file exists
    [dataFile, ~, dataFileExists] = coastfilename(mfilename, ...
        'Upscale', upscale, 'LatLim', latlim, ...
        'Buffer', buf, 'MoreBuffers', moreBuf);

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

    % Find the coastline
    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBuf, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    % Subtract the coast and earthquake mask from the ocean
    oceanPoly = subtract(oceanPoly, coastPoly);

    XY = closecoastline(oceanPoly.Vertices);

    % Think twice before upscaling!
    % Bezier splines may lead to holes between oceans
    if upscale ~= 0 && upscale ~= 1
        XY = bezier(XY, upscale);
    end

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
        'Latlim', [-16, -15], 'Lonlim', [127, 129], ...
        'LongitudeOrigin', lonOrigin);

    if buf >= 0.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [13, 28; 24, 30; 10, 13; 1, 3.5], ...
            'Lonlim', [34, 45; 45, 56; 43, 45.5; 100, 103], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [10, 14; -9.5, 8; 4, 8], ...
            'Lonlim', [45, 48; 115, 116; 98, 100], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 1.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [6, 8], 'Lonlim', [97, 99], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 3.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-24, -20], 'Lonlim', [34, 40], ...
            'LongitudeOrigin', lonOrigin);
    end

end
