%% OCEANS
% Finds the coordinates of some of the worlds' oceans so we can combine
% its localization kernel with those for the missing continents to turn
% into Slepian eigenfunctions for all of the world's oceans.
%
% Syntax
%   XY = oceans(Upscale, Buffer)
%   oceans(__)
%
% Input arguments
%   Upscale - How many times to upscale the data
%       This option should be used carefully, as spline interpolation do not work well at sharp corners, and the original sampling rate should be sufficient for most applications.
%       The default value is 0 (no upscaling).
%   Buffer - The (negative) buffer from the coastlines in degrees
%       The buffer is specified as a scalar or a 1-by-1+2n cell.
%       The default value is 0.
%   Inclang - The latitudes of the polar caps in degrees
%       The inclination angle can be specified as a scalar or a 1-by-2 vector.
%       - If a scalar, the inclination angle is applied to both the north and south polar caps.
%       - If a 1-by-2 vector, the first element is the inclination angle of the north polar cap, and the second element is the inclination angle of the south polar cap.
%       The default value is 90 (no polar caps).
%   ForceReload - Force the function to reload the data
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the Atlantic ocean
%
% Last modified by williameclee-at-arizona.edu, Jul 2nd, 2024

function varargout = oceans(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, inclang, buf, moreBuf, forceReload, lonOrigin] = ...
        parsecoastinputs(varargin, 'DefaultLongitudeOrigin', 180);
    oceanTitle = 'All Ocean';
    oceanParts = ...
        {'Atlantic Ocean', ...
         'Indian Ocean', ...
         'Pacific Ocean, eastern part', ...
     'Pacific Ocean, western part'};

    %% Check if the data file exists
    [dataFile, ~, dataFileExists] = coastfilename(mfilename, ...
        'Upscale', upscale, 'Inclang', inclang, ...
        'Buffer', buf, 'MoreBuffer', moreBuf);

    if dataFileExists && ~forceReload
        load(dataFile, 'XY')

        varargout = returncoastoutput(nargout, XY, oceanTitle);

        return
    end

    %% Compute the data
    [oceanPoly, oceanLatlim, oceanLonlim] = ...
        findoceanboundary(oceanParts, inclang, lonOrigin);

    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBuf, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    oceanPoly = subtract(oceanPoly, coastPoly);

    XY = closecoastline(oceanPoly.Vertices);

    % Convert the data to the right format
    [X, Y] = poly2cw(XY(:, 1), XY(:, 2));
    % [Y, X] = flatearthpoly(Y, X, lonOrigin);
    XY = removeduplicatevertices([X, Y]);

    %% Save and return required data
    save(dataFile, 'XY', '-v7.3')
    fprintf('%s saving %s\n', upper(mfilename), dataFile)

    varargout = returncoastoutput(nargout, XY, oceanTitle);

end

function coastPoly = manualadjustment(coastPoly, buf, lonOrigin)

    coastPoly = addlandregion(coastPoly, ...
        'Latlim', [-16, -15], 'Lonlim', [127, 129], ...
        'LongitudeOrigin', lonOrigin);

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-40, -39], 'Lonlim', [143.5, 145.5], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-40, -39.9], 'Lonlim', [147, 147.2], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [18, 19], 'Lonlim', [-76, -75], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [45, 50], 'Lonlim', [-65, -60], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [58, 61], 'Lonlim', [156, 160], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [38, 39], 'Lonlim', [122, 124], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [22, 24], 'Lonlim', [118, 120], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [9, 11], 'Lonlim', [124, 125], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [14, 24], 'Lonlim', [35, 45], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [24, 30], 'Lonlim', [45, 55], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [11, 12], 'Lonlim', [44, 45.5], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 2
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [37, 43], 'Lonlim', [130, 138], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [54, 58], 'Lonlim', [0, 6], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 3
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [20, 30], 'Lonlim', [-95, -85], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [10, 20], 'Lonlim', [-85, -70], ...
            'LongitudeOrigin', lonOrigin);
    end

end
